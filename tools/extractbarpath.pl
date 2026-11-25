#!/usr/bin/perl
# Extract strains/stresses along a chain of lattice3d (or old latticeframe3d)
# bar elements in OOFEM, fast.
#
# For each step i = 1..nSteps, reads <gp_basename>.<i>.gp (gpexportModule,
# vars 138 139 ...) and writes <output_prefix><i>.dat with columns
#   coord (arc length of element midpoint along the chain), strain, force
# for the elements on the shortest bar path from start_node to end_node.
use strict;
use warnings;

my $usage = <<"END";
Usage: $0 <gp_basename> <oofem_input> <start_node> <end_node> <n_steps> <output_prefix>
  gp_basename    e.g. oofem.out.m1  (files oofem.out.m1.<i>.gp)
  oofem_input    OOFEM input file with the mesh, e.g. oofem.in
  start_node     node number at one end of the bar chain
  end_node       node number at the other end of the bar chain
  n_steps        number of steps to process (1..n_steps)
  output_prefix  output files are <output_prefix><i>.dat
END

die "Error: expected 6 arguments, got " . scalar(@ARGV) . ".\n$usage"
  unless @ARGV == 6;
my ($gp_base, $meshFile, $startNode, $endNode, $nSteps, $outputName) = @ARGV;

for ([start_node => $startNode], [end_node => $endNode], [n_steps => $nSteps]) {
    my ($name, $val) = @$_;
    die "Error: $name must be a positive integer, got '$val'.\n$usage"
      unless $val =~ /^\d+$/ && $val > 0;
}

# ----------- READ MESH (ONCE) -----------
open my $IN, "<", $meshFile or die "Could not open $meshFile\n";

my ($nNodes, $nElems);
my (@nodes_x, @nodes_y, @nodes_z);      # dense idx 0..nNodes-1
my %node_dense;                          # OOFEM nodeNumber -> dense idx

# bar element storage (lattice3d or latticeframe3d)
my @lf_n1; my @lf_n2; my @lf_num;        # dense latticeframe elements
my %lf_by_num;                           # elementNumber -> dense idx (bars only)
my %adj;                                 # nodeNumber -> [neighborNode, lf_eid]
my %edge_by_pair;                        # "min,max" -> lf_eid (to map node pair to element)

while (my $line = <$IN>) {
    if ($line =~ /ndofman/i) {
        my @t = split /\s+/, $line;
        for my $i (0..$#t-1) {
            $nNodes = $t[$i+1] if $t[$i] =~ /ndofman/i;
            $nElems = $t[$i+1] if $t[$i] =~ /nelem/i;
        }
        die "Could not read ndofman/nelem in $meshFile\n" unless defined $nNodes && defined $nElems;

        # ---- nodes block ----
        for (my $k = 0; $k < $nNodes; ++$k) {
            my $ln = <$IN>;
            my @f  = split /\s+/, $ln;
            # Original script used positions [1]=nodeNumber, [4..6]=x,y,z
            my $nodeNumber = $f[1];
            my ($x,$y,$z)  = @f[4,5,6];
            my $dense = scalar @nodes_x;
            $node_dense{$nodeNumber} = $dense;
            $nodes_x[$dense] = $x;
            $nodes_y[$dense] = $y;
            $nodes_z[$dense] = $z;
        }

        # ---- elements block (filter to bar elements) ----
        for (my $k = 0; $k < $nElems; ++$k) {
            my $ln = <$IN>;
            my @f  = split /\s+/, $ln;
            next unless @f; # skip blank lines
            # Require lattice3d (or old latticeframe3d) as element type:
            next unless $f[0] =~ /^lattice(frame)?3d$/i;

            # Keep same field positions as your original:
            my $elnum = $f[1];
            my ($n1,$n2) = @f[4,5];

            my $eid = scalar @lf_n1;
            $lf_n1[$eid]  = $n1;
            $lf_n2[$eid]  = $n2;
            $lf_num[$eid] = $elnum;
            $lf_by_num{$elnum} = $eid;

            push @{ $adj{$n1} }, [$n2, $eid];
            push @{ $adj{$n2} }, [$n1, $eid];

            my ($a,$b) = ($n1<$n2) ? ($n1,$n2) : ($n2,$n1);
            $edge_by_pair{"$a,$b"} = $eid;
        }
        last; # done after first mesh block
    }
}
close $IN;

die "Start node $startNode not present in bar (lattice3d) graph.\n"
  unless exists $adj{$startNode};
die "End node $endNode not present in bar (lattice3d) graph.\n"
  unless exists $adj{$endNode};

# ----------- BFS PATH OVER BAR ELEMENT GRAPH -----------
my %visited; my %parent;
my @q = ($startNode);
$visited{$startNode} = 1;
PUSH: while (@q) {
    my $u = shift @q;
    last if $u == $endNode;
    for my $p (@{ $adj{$u} || [] }) {
        my ($v, undef) = @$p;
        next if $visited{$v};
        $visited{$v} = 1;
        $parent{$v}  = $u;
        push @q, $v;
    }
}
die "No path between $startNode and $endNode using bar (lattice3d) elements.\n"
  unless $visited{$endNode};

# Reconstruct node path
my @path_nodes;
for (my $v = $endNode; defined $v; $v = $parent{$v}) {
    unshift @path_nodes, $v;
    last if $v == $startNode;
}

# Convert node path to element indices & mid-arc coordinates
my @path_eid;            # dense bar element indices
my @path_elnum;          # element numbers (as in .gp)
my @path_midcoord;       # s(midpoint) along chain
{
    my $s = 0.0;
    for (my $i=0; $i<$#path_nodes; ++$i) {
        my ($a,$b) = ($path_nodes[$i], $path_nodes[$i+1]);
        my ($m,$n) = ($a<$b) ? ($a,$b) : ($b,$a);
        my $eid = $edge_by_pair{"$m,$n"};
        die "Internal error mapping edge $a-$b to element.\n" unless defined $eid;
        push @path_eid,   $eid;
        push @path_elnum, $lf_num[$eid];

        my $i1 = $node_dense{$a};
        my $i2 = $node_dense{$b};
        my $dx = $nodes_x[$i1] - $nodes_x[$i2];
        my $dy = $nodes_y[$i1] - $nodes_y[$i2];
        my $dz = $nodes_z[$i1] - $nodes_z[$i2];
        my $L  = sqrt($dx*$dx + $dy*$dy + $dz*$dz);

        push @path_midcoord, $s + 0.5*$L;
        $s += $L;
    }
}

# Membership set for quick filtering when reading .gp
my %want = map { $_ => 1 } @path_elnum;

# ----------- PROCESS GP FILES -----------
for my $i (1..$nSteps) {
    my $gp = "$gp_base.$i.gp";

    unless (-e $gp) {
        warn "Warning: file $gp not found, skipping step $i\n";
        next;
    }

    open my $GP, "<", $gp or do { 
        warn "Could not open $gp, skipping\n";
        next;
    };

    # Skip the header like your original
    my $h1 = <$GP>;
    my $h2 = <$GP>;
    my $h3 = <$GP>;
    my $h4 = <$GP>;

    my %vals; # elnum -> [strain, stress]
    while (my $ln = <$GP>) {
        next if $ln =~ /^\s*$/;
        my @f = split /\s+/, $ln;
        next unless @f >= 15;
        next unless $f[0] =~ /^\d+$/;
        my $el = $f[0];
        next unless $want{$el};
        my $strain = $f[10];
        my $stress = $f[14];
        $vals{$el} = [$strain, $stress];
    }
    close $GP;

    my $out = "$outputName$i.dat";
    open my $OUT, ">", $out or do {
        warn "Cannot write $out\n";
        next;
    };
    print $OUT "# coord strain stress\n";
    for my $k (0..$#path_elnum) {
        my $el = $path_elnum[$k];
        next unless exists $vals{$el};
        printf $OUT ("%g %s %s\n", $path_midcoord[$k], $vals{$el}[0], $vals{$el}[1]);
    }
    close $OUT;
}
