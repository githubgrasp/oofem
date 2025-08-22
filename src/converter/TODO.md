TODO — Converter & Generator

1) Output refactor (modular writers)
	•	Create src/converter/output/ with one file per case:
	•	Writer.hpp (pure-virtual base)
	•	WriterFactory.hpp/.cpp (string → writer)
	•	Writer3DSMTM.cpp, Writer3DTM.cpp, Writer3DSM.cpp, … (one per gridType)
	•	Move each giveXXXOutput(...) body from Grid into a dedicated Writer* class.
	•	In Grid::giveOutput(name): auto w = WriterFactory::make(gridType); w->write(*this, name);
	•	Keep legacy filenames/format intact (byte-for-byte if possible).

2) Standalone test harness (separate from OOFEM tests)
	•	Add tests/ with structure:
	•	tests/data/ (tiny inputs: minimal, periodic, cylinder, edge-cases)
	•	tests/generator/ & tests/converter/ (test code)
	•	Use CTest (already in CMake):
	•	Add option(USE_GEN_TESTS ON) / option(USE_CONV_TESTS ON)
	•	add_test(NAME gen_minimal COMMAND generator.exec tests/data/minimal.in)
	•	Post-run scripts to verify:
	•	file exists (e.g., oofem.in, .vtu)
	•	line-count / regex sanity checks (no “undefined id”, etc.)
	•	Provide a tiny golden-output compare:
	•	tests/scripts/diff_trim.py (ignores whitespace, float tol via regex)

3) CI & build quality
	•	GitHub Actions workflow:
	•	matrix: {macOS, ubuntu} (release+asan)
	•	cache Homebrew/apt deps (qhull)
	•	run CTest and upload artifacts on failure
	•	Add ASan/UBSan config:
	•	-fsanitize=address,undefined -fno-omit-frame-pointer -g
	•	CMake option: -DUSE_ASAN=ON

4) Input parsing robustness
	•	Case-insensitive keywords (strcasecmp path):
	•	Normalize to lowercase at parse boundary.
	•	Ensure IntArray::resizeWithValues used where appending is intended.
	•	Replace remaining raw fopen/printf with converter::fopen_or_die & std::string.

5) Performance checkpoints
	•	Add timing macros around: point generation, octree insert, neighbor search, file write.
	•	Optional: skip mirroring for far-boundary candidates (config flag).

6) Code health & safety
	•	Ensure every new-allocated object is owned by a single container and freed in Grid dtor (already mostly done).
	•	Verify GridLocalizer::init(...) is called exactly once and after bbox is valid.
	•	Replace any growTo remnants with resize (or helpers ensure_size1).
	•	Avoid shadowed virtuals: match base signatures (const-correctness, params).

7) Docs (lightweight for now)
	•	Create docs/ with:
	•	README.md (build/run, minimal examples)
	•	inputs.md (keyword list, case-insensitive note, examples)
	•	dev.md (layout, how Writers work, test rules)

8) Developer helpers
	•	scripts/format.sh (uncrustify / clang-format with updated options)
	•	scripts/run_gen.sh, scripts/run_conv.sh (with caffeinate on macOS)
	•	scripts/profile.sh (simple perf/instruments/time wrapper)