# XML input format

The XML input format can be used everywhere where the text-based format is used. Currently, it requires oofem to be compiled with `-DUSE_XML=1`, which will links against system-installed [pugixml](https://pugixml.org) library for parsing XML. XML input will be activated for input files having the `.xml` extension automatically.

- records (lines in the text input) are expressed as entities, fields are entity attributes;
- entities are nested;
- text between elements (`PCDATA` as it is called in XML) is only used for `<Output>...</Output>` and `<Description>...</Description>` (not anywhere else);
- there is an [experimental tool](https://github.com/eudoxos/oofem-xml/) to convert from text to XML which currently passes (almost) all unit tests; that repository also contains all tests converted to XML;
- tags and attributes are case-insensitive, but message (warning) will be shown for non-matching case;
- special characters in field names are replaced by `_` to conform with XML attribute names (thus `f(t)` becomes `f_t_`, `a/c` becomes `a_c`)

In comparison to the text format:

- XML is machine readable and writeable, using usual XML tools;
- versioned (`<oofem version="..."`), with automated upgrade paths;
- number of enumerated sub-itmes (such as `nelem` in a domain) are inferred (the {ref}`ComponentsSizeRecord` is therefore completely missing);
- ordering of heterogeneous records (such as whether `<Sets>` are specified before `<Elements>`) is arbitrary;
- `id`'s of enumerated elements are optional (but checked);
- XML comments `<!-- .... -->` are supported;
- XML fragments can be included via `<xi:include>` directive;
- `#%BEGIN_CHECK% ... #%END_CHECK%` (`ErrorCheckingExportModule`) is specified as XML, using `<ErrorCheck>` group, instead of special syntax.


## Field formats

Since attributes in XML are string-only, here we describe how other types are represented:

:::{list-table}
:header-rows: 1

- - type
  - description
  - notes
- - integers, floats
  - usual ASCII representation
  -
- - flags
  - attributes without value are written as `attribute=""` in XML
  - non-empty value will produce warning about not having been processed
- - boolean
  - true (`1`, `y`, `Y` and a few others), false (`0`, `n`, `N` and others)
  - unrecognized value produces error
- - 1D arrays (ints, floats)
  - space-separated
  -
- - 2D matrices
  - colon-separated rows, space-separated elements
  - e.g. `a11 a21; a12 a22`
- - dictionaries
  - colon-separated items, space-separated key and value
  - e.g `key1 value; key2 value`
- - ranges
  - comma-separated items or contiguous ranges with dash (minus)
  - e.g. `1-4,7,9-11`; as alternative: `1..4,7,9..1`
- - enums
  - specified as number of keyword
  -
:::


## Input example


```{code-block} xml
<!-- the top-level tag must be <oofem ...> -->
<oofem version="1">
  <Output>spring01.out</Output>
  <Description>Patch test of 1D spring element along x axis</Description>
  <Analysis type="StaticStructural" nsteps="1" />
  <Domain domain="1dtruss">
    <!-- flags have empty value -->
    <OutputManager tstep_all="" dofman_all="" element_all="" />
    <!-- id params are optional here, things are numbered from 1 by default -->
    <Nodes>
      <Node id="1" coords="0. 0. 0."  />
      <Node id="2" coords="0. 0. 0."/>
    </Nodes>
    <Elements>
      <Spring nodes="2 1 2" crosssect="1" mode="0" k="2.0"/>
    </Elements>
    <CrossSections>
      <SimpleCS thick="0.1" width="1.0" material="1" set="1"/>
    </CrossSections>
    <Materials>
      <dummymat />
    </Materials>
    <BoundaryConditions>
      <BoundaryCondition id="1" loadtimefunction="1" dofs="1" values="0.0" set="2"/>
      <NodalLoad id="2" loadtimefunction="1" dofs="1" components="1.0" set="3"/>
    </BoundaryConditions>
    <InitialConditions/>
    <LoadTimeFunctions>
      <ConstantFunction f_t_="1.0"/>
    </LoadTimeFunctions>
    <Sets>
      <Set nodes="1"/>
      <Set nodes="1"/>
      <Set id="3" nodes="2"/>
    </Sets>
  </Domain>
  <ExportModules>
    <errorcheck tolerance="1e-4">
      <NODE tStep="1" number="2" dof="1" unknown="d" value="0.5"/>
      <REACTION tStep="1" number="1" dof="1" value="-1"/>
    </errorcheck>
  </ExportModules>
</oofem>
```
