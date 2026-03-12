# lib/matrix — AI Context

## Purpose

`matrix` is Layer 2 in Ka-Boost. It provides two distinct surfaces: (1) procedural 1D array operations on `ARRAY[*] OF REAL` (init, arithmetic, linspace, shift), and (2) a GPP-expanded 2D matrix class template (`2dArray.klc`) that generates named Karel programs for arbitrarily-sized matrices with full linear algebra (multiplication, inverse, determinant, cofactor, trace). Karel has no built-in matrix or array algebra; this module fills that gap. The 4×4 variant is the backbone of `lib/pose/matpose`, which implements homogeneous coordinate transforms for the 5-axis slicer pipeline.

## Repository Layout

```
lib/matrix/
├── package.json                         -- rossum manifest (depends: math, errors, ktransw-macros)
├── readme.md                            -- sparse developer note
├── src/
│   └── matrix.kl                        -- compiled Karel program; all matrix__ routines live here
├── include/
│   ├── matrix.klh                       -- public declarations for matrix__ functions + short aliases
│   ├── matrix.klt                       -- OP_* preprocessor constants (OP_ZEROS..OP_SMULT)
│   ├── matrix.private.klh               -- private declarations (operations_1D, init_1D)
│   └── classes/
│       ├── 2dArray.klh                  -- class header (declare_member for all public methods)
│       ├── 2dArray.klc                  -- class implementation (GPP template; expands to named program)
│       ├── 2dArray.private.klh          -- private method declarations (init_, operations_, dot_, etc.)
│       └── templates/
│           ├── carr3.klt                -- 3×3 REAL (rotation matrices)
│           ├── carr4.klt                -- 4×4 REAL (homogeneous transforms)
│           ├── carr10.klt               -- 10×10 REAL (general)
│           ├── carr23.klt               -- 2×3 INTEGER (general)
│           └── carr305.klt              -- 30×5 REAL (path data)
└── test/
    ├── test_matrix.kl                   -- KUnit tests for 1D functions
    └── test_matrix2.kl                  -- KUnit tests for all five 2D class variants
```

## Full API Reference

### 1D Operations (`src/matrix.kl`, declared in `include/matrix.klh`)

All 1D routines operate on `ARRAY[*] OF REAL`. Short aliases (via `declare_function`) are listed after each signature.

| Routine | Alias | Signature | Behavior |
|---------|-------|-----------|----------|
| `matrix__full_1D` | `mtrx__full1R` | `(arr : ARRAY[*] OF REAL; number : REAL)` | Fill every element with `number` |
| `matrix__zeros_1D` | `mtrx__zros1R` | `(arr : ARRAY[*] OF REAL)` | Fill every element with `0.0` |
| `matrix__ones_1D` | `mtrx__ones1R` | `(arr : ARRAY[*] OF REAL)` | Fill every element with `1.0` |
| `matrix__random_1D` | `mtrx__rand1R` | `(arr : ARRAY[*] OF REAL)` | Fill with `math__rand` [0,1]; seeds with `GET_TIME * i` each element |
| `matrix__linspace` | `mtrx__lspce_r` | `(start_idx, stop_idx : REAL; arr : ARRAY[*] OF REAL; endpoint : BOOLEAN)` | Evenly-spaced values. `endpoint=TRUE` → includes `stop_idx`; `endpoint=FALSE` → excludes it (NumPy-style) |
| `matrix__eye` | `mtrx__eye` | `(arr : ARRAY[*,*] OF REAL; row : ARRAY[*] OF REAL)` | Write identity into 2D native array. `row` is a dummy 1D array sized to the column count — Karel cannot query second-dimension length of `ARRAY[*,*]` |
| `matrix__shift_s` | `mtrx__shifts` | `(arr : ARRAY[*] OF STRING; amount : INTEGER; updown : BOOLEAN; unint_arr : ARRAY[*] OF STRING)` | Cyclic rotate STRING array by `amount`. `updown=TRUE` → shift up (higher indices); `FALSE` → shift down. `unint_arr` is required temp buffer same size as `arr`; result is written back to `arr` |
| `matrix__shift_r` | `mtrx__shiftr` | `(arr : ARRAY[*] OF REAL; amount : INTEGER; updown : BOOLEAN; unint_arr : ARRAY[*] OF REAL)` | Same as `shift_s` but for REAL arrays |
| `matrix__add_1D` | `mtrx__add1R` | `(arr1, arr2 : ARRAY[*] OF REAL; out_arr : ARRAY[*] OF REAL)` | Element-wise `arr1 + arr2`. Arrays must be same length; UNINIT elements are passed through from the initialized operand |
| `matrix__sub_1D` | `mtrx__sub1R` | `(arr1, arr2 : ARRAY[*] OF REAL; out_arr : ARRAY[*] OF REAL)` | Element-wise `arr1 - arr2` |
| `matrix__mul_1D` | `mtrx__mul1R` | `(arr1, arr2 : ARRAY[*] OF REAL; out_arr : ARRAY[*] OF REAL)` | Element-wise `arr1 * arr2` |
| `matrix__div_1D` | `mtrx__div1R` | `(arr1, arr2 : ARRAY[*] OF REAL; out_arr : ARRAY[*] OF REAL)` | Element-wise `arr1 / arr2`. Division by zero is silently skipped: `out_arr[i] = arr1[i]` when `arr2[i] = 0` |

### OP_* Constants (`include/matrix.klt`)

Preprocessor defines used internally by dispatch routines. Not part of the public API, but visible to any file that `%include matrix.klt`.

```
OP_ZEROS    = 1
OP_EYE      = 2
OP_FULL     = 3
OP_RANDOM   = 4
OP_ONES     = 5
OP_ADD      = 6
OP_SUBTRACT = 7
OP_MULTIPLY = 8   -- 1D element-wise multiply
OP_DIVIDE   = 9
OP_SMULT    = 8   -- 2D scalar multiply (same integer as OP_MULTIPLY — intentional, different dispatch paths)
```

### 2D Array Class (`include/classes/2dArray.klc`)

The class template is instantiated by expanding `2dArray.klc` via GPP with a `.klt` config that defines:
- `type_name` — name for the generated 2D array type (e.g. `t_carr4`)
- `t_rows`, `t_columns` — integer dimensions
- `inner_type` — `REAL` or `INTEGER`
- `UNIT_TESTING` — define to include `test_row` / `testing` routines

After expansion, all method names follow `class_name__method` (e.g. `carr4__mult`). Short 4-character aliases are also generated (e.g. `carr4__mult`).

**Type macros generated by `define_type.m`:**
- `arry(type_name)` — the 2D row-of-rows type (e.g. `ARRAY[4] OF t_carr4`)
- `arrx(type_name)` — the 1D row type (e.g. `ARRAY[4] OF REAL` for a 4-column matrix)

#### Creation / Parsing

| Method | Alias | Signature | Behavior |
|--------|-------|-----------|----------|
| `create_array_from_string` | `cafs` | `(sarr : STRING; col_delim : STRING; row_delim : STRING) : arry(type_name)` | Parse delimited string into matrix. Splits on `row_delim` first, then `col_delim`. Type-converts via `s_to_rarr` (REAL) or `s_to_iarr` (INTEGER) |
| `create_row_from_string` | `crws` | `(sarr : STRING; col_delim : STRING) : arrx(type_name)` | Parse one row from a delimited string |

#### Clearing / Setting

| Method | Alias | Signature | Behavior |
|--------|-------|-----------|----------|
| `clear` | `clr` | `(arr : arry(type_name))` | Set all elements to UNINIT via uninitialized local variable |
| `set_row` | `setr` | `(row_idx : INTEGER; row : arrx(type_name); out_arr : arry(type_name))` | Replace row at index. Aborts if `row_idx > ARRAY_LEN(out_arr)` |
| `set_from_array` | `star` | `(in_arr : ARRAY[*,*] OF inner_type; out_arr : arry(type_name))` | Copy from native 2D array to class array. Arrays must be same row count |
| `set_row_from_array` | `srfa` | `(in_row : ARRAY[*] OF inner_type; row_idx : INTEGER; out_arr : arry(type_name))` | Copy native 1D array into a specific row |

#### Getting

| Method | Alias | Signature | Notes |
|--------|-------|-----------|-------|
| `get_col` | `getc` | `(col_idx : INTEGER; arr : arry(type_name)) : arrx(type_name)` | **Square matrices only** (guarded by `%ifeq t_rows t_columns`). Transposes internally to extract column |

#### Conversion

| Method | Alias | Signature | Behavior |
|--------|-------|-----------|----------|
| `convert_to_array` | `cvar` | `(in_arr : arry(type_name); out_arr : ARRAY[*,*] OF inner_type)` | Copy class array to native 2D Karel array. **Note:** Both dimension checks compare `ARRAY_LEN(out_arr)` (first dimension only) — this is a bug for non-square matrices |

#### Initialization

| Method | Alias | Signature | Behavior |
|--------|-------|-----------|----------|
| `full` | `full` | `(number : inner_type; out_arr : arry(type_name))` | Fill all elements with `number` |
| `zeros` | `zero` | `(out_arr : arry(type_name))` | Fill with 0 |
| `ones` | `ones` | `(out_arr : arry(type_name))` | Fill with 1 |
| `eye` | `eye` | `(out_arr : arry(type_name))` | Identity matrix (1 on diagonal, 0 elsewhere) |
| `random` | `rand` | `(out_arr : arry(type_name))` | Fill with random values; INTEGER uses `math__rand_int`, REAL uses `math__rand` [0,1]. Seeds each cell with `GET_TIME * i * j` |

#### Element-wise and Scalar Operations

| Method | Alias | Signature | Behavior |
|--------|-------|-----------|----------|
| `add` | `add` | `(arr1, arr2 : arry(type_name)) : arry(type_name)` | Element-wise addition. UNINIT elements fall back to the initialized operand |
| `subtract` | `subt` | `(arr1, arr2 : arry(type_name)) : arry(type_name)` | Element-wise subtraction |
| `smult` | `smlt` | `(arr : arry(type_name); val : inner_type) : arry(type_name)` | Scalar multiplication — all elements × `val` |

#### Square Matrix Operations (only compiled when `t_rows = t_columns`)

| Method | Alias | Signature | Behavior |
|--------|-------|-----------|----------|
| `mult` | `mult` | `(arr1, arr2 : arry(type_name)) : arry(type_name)` | Matrix multiplication. Transposes `arr2` internally, then dot-products rows |
| `transpose` | `trnp` | `(arr : arry(type_name)) : arry(type_name)` | Returns transposed copy |
| `inverse` | `inv` | `(arr : arry(type_name)) : arry(type_name)` | Computes `(1/det) * transpose(cofactor(arr))`. Aborts if `abs(det) < 0.05` |
| `det` | `det` | `(arr : arry(type_name)) : inner_type` | Recursive cofactor expansion along first row. Base case: 1×1 returns `arr[1,1]` |
| `cofactor` | `cfct` | `(arr : arry(type_name)) : arry(type_name)` | Full cofactor matrix. Each element: `(-1)^(i+j) * det(minor)` |
| `trace` | `trac` | `(arr : arry(type_name)) : inner_type` | Sum of diagonal elements |

#### Testing (only compiled when `UNIT_TESTING` is defined)

| Method | Alias | Signature | Behavior |
|--------|-------|-----------|----------|
| `test_row` | `trw` | `(expected, actual : arrx(type_name)) : BOOLEAN` | Fuzzy-compare one row; EPSILON = 0.005. Writes error message to KUnit global |
| `testing` | `test` | `(expected, actual : arry(type_name)) : BOOLEAN` | Compare entire matrix row-by-row using `test_row` |

### Built-in Template Configs

| File | `type_name` | Dims | `inner_type` | Use |
|------|------------|------|-------------|-----|
| `carr3.klt` | `t_carr3` | 3×3 | REAL | Rotation matrices |
| `carr4.klt` | `t_carr4` | 4×4 | REAL | Homogeneous transforms |
| `carr10.klt` | `t_carr10` | 10×10 | REAL | General |
| `carr23.klt` | `t_carr23` | 2×3 | INTEGER | General |
| `carr305.klt` | `t_carr305` | 30×5 | REAL | Path data |

## Core Patterns

### Pattern 1: Instantiating a 2D Matrix Class

The `%class` directive expands `2dArray.klc` into a new Karel program named `class_name`. The `.klt` config must be included before the expansion so GPP can substitute the macros.

```
-- In your .kl file:
%include define_type.m
%include carr4.klt
t_arr2d_ref(t_rows, t_columns, type_name, inner_type)   -- declares t_carr4 type

%class carr4('2dArray.klc','2dArray.klh','carr4.klt')   -- generates carr4 program
```

After this, all methods are accessible as `carr4__zeros(mat)`, `carr4__mult(a, b)`, etc. The program `carr4.pc` is compiled and deployed alongside your program.

### Pattern 2: Custom Matrix Type for a Module

When the built-in templates don't fit, create a custom `.klt`:

```
-- mymod/mymat.klt
%include define_type.m
%defeval type_name t_mymat
%defeval t_rows 4
%defeval t_columns 4
%defeval inner_type REAL
-- omit UNIT_TESTING for production builds
```

Then in your `.kl`:
```
%include mymat.klt
t_arr2d_ref(t_rows, t_columns, type_name, inner_type)
%class cmat('2dArray.klc','2dArray.klh','mymat.klt')
```

This is the pattern used by `lib/pose/matpose` (`matarr.klt` → `cmat` program).

### Pattern 3: Building a Rotation Matrix with `set_row` (from `lib/pose/matpose/matpose.kl`)

```
ROUTINE matpose__rotx
  VAR
    row : arrx(type_name)
    out_mat : arry(type_name)
  BEGIN
    row[1] = 1 ; row[2] = 0 ; row[3] = 0 ; row[4] = 0
    cmat__set_row(1, row, out_mat)
    row[1] = 0 ; row[2] = COS(angle) ; row[3] = -1*SIN(angle) ; row[4] = 0
    cmat__set_row(2, row, out_mat)
    row[1] = 0 ; row[2] = SIN(angle) ; row[3] = COS(angle) ; row[4] = 0
    cmat__set_row(3, row, out_mat)
    row[1] = 0 ; row[2] = 0 ; row[3] = 0 ; row[4] = 1
    cmat__set_row(4, row, out_mat)
    RETURN(out_mat)
  END matpose__rotx
```

### Pattern 4: Parsing a Matrix from a String

Useful for loading config or CSV data into a matrix:

```
-- from test_matrix2.kl:
VAR
  mat : arry(type_name)
  row : arrx(type_name)
BEGIN
  mat = carr3__create_array_from_string('2,-3,9;2,0,-1;1,4,5', ',', ';')
  -- mat[1] = [2, -3, 9], mat[2] = [2, 0, -1], mat[3] = [1, 4, 5]

  row = carr23__create_row_from_string('1,2,4', ',')
  -- row = [1, 2, 4]  (INTEGER, from carr23 config)
```

### Pattern 5: 1D linspace for Evenly-Spaced Sampling

```
-- from test_matrix.kl:
VAR
  pts : ARRAY[5] OF REAL
BEGIN
  -- endpoint=TRUE: includes stop value
  matrix__linspace(2.0, 3.0, pts, TRUE)
  -- pts = [2.0, 2.25, 2.5, 2.75, 3.0]

  -- endpoint=FALSE: excludes stop value (NumPy default)
  matrix__linspace(2.0, 3.0, pts, FALSE)
  -- pts = [2.0, 2.2, 2.4, 2.6, 2.8]
```

## Common Mistakes

| Mistake | Symptom | Fix |
|---------|---------|-----|
| Calling `mult`, `det`, `inverse`, `transpose`, `get_col`, `trace`, `cofactor` on a non-square class | Compile error: routine not found | These are compiled only when `t_rows = t_columns`. Create a square template or restructure the algorithm |
| Forgetting the `unint_arr` temp buffer in `matrix__shift_s` / `matrix__shift_r` | Runtime error or corrupted data | Declare a same-sized array of the same type and pass it as the 4th argument; it is modified in-place |
| Using `matrix__eye` on a 2D class array (instead of a native `ARRAY[*,*]`) | Type mismatch compile error | `matrix__eye` is the 1D procedural version for native `ARRAY[*,*]`; use `carr4__eye(mat)` for class arrays |
| Setting `inner_type = STRING` in a custom `.klt` | Compile error at GPP `%error` directive inside `2dArray.klc` | Only `REAL` and `INTEGER` are supported as `inner_type` |
| Calling `inverse` on a near-singular matrix | ABORT error: "Matrix is not invertable" | The singularity threshold is `abs(det) < 0.05` — much larger than floating-point epsilon. Well-conditioned matrices only |
| Forgetting `t_arr2d_ref(...)` after including the `.klt` | Type `t_carr*` undefined; compile error | `t_arr2d_ref` expands the `ARRAY` type declaration from `define_type.m`; it must appear in the `VAR` section or top of the program |
| `OP_SMULT = 8 = OP_MULTIPLY` confusion | Using `OP_SMULT` in a 1D context hits the `ELSE` branch → ABORT | Do not call `matrix__operations_1D` directly with `OP_SMULT`; use the public `matrix__mul_1D` for element-wise 1D multiply |
| `convert_to_array` on non-square matrices | Spurious `ARR_LEN_MISMATCH` abort | The second `ARRAY_LEN(out_arr)` check compares first-dimension length against `t_columns`; it only works correctly for square matrices |

## Dependencies

### This module depends on:
- `math` — `math__srand`, `math__rand`, `math__rand_int`
- `errors` — `karelError`, `CHK_STAT`, `ARR_LEN_MISMATCH`, `INVALID_INDEX`, `VAL_OUT_OF_RNG`
- `ktransw-macros` — `declare_function`, `declare_member`, `namespace.m`, `header_guard.m`, `define_type.m`
- `Strings` — `s_to_arr`, `s_to_iarr`, `s_to_rarr`, `i_to_s`, `r_to_s` (used inside `2dArray.klc`)
- `KUnit` — test routines only (`kunit.klh`)

### Modules that depend on this module:
- `lib/pose` — sole external consumer. `matpose` uses `2dArray.klc` with custom `matarr.klt` (4×4) and `rotarr.klt` (3×3). `pose.kl` imports `matrix.klh` for 1D routines. All homogeneous transform math in the 5-axis pipeline flows through `cmat` and `crot` class instances.

## Build / Integration Notes

- `matrix.kl` compiles to `matrix.pc` (the 1D routine program). Every 2D class variant (`carr3.pc`, `carr4.pc`, etc.) compiles to its own `.pc` file because each `%class` expansion generates a separate Karel program.
- When using `2dArray.klc` in a module, include `matrix` as a dependency in `package.json`. rossum will resolve it and add `matrix.pc` to the build.
- The `UNIT_TESTING` define gates `test_row` and `testing` methods. These are enabled in all bundled templates (`carr*.klt`) but not in pose's `matarr.klt` (production config has no `%define UNIT_TESTING`). Strip `UNIT_TESTING` in production `.klt` files to avoid dead code on the controller.
- The `2dArray.klc` `%include kunit.klh` is unconditional, meaning `KUnit` is always a transitive dependency even in production builds. This is a minor build hygiene issue but does not affect runtime behavior if `kunit.pc` is not running.
- `get_col`, `transpose`, `mult`, `det`, `cofactor`, `inverse`, `trace` are conditionally compiled with `%ifeq t_rows t_columns`. GPP evaluates this as a string comparison — `t_rows` and `t_columns` must expand to the same decimal integer string, not equivalent expressions.
