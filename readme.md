# matrix

FANUC Karel 1D array operations and GPP-templated 2D matrix class for linear algebra in the Ka-Boost toolchain.

## Overview

Karel has no native array algebra and no generics. `matrix` fills both gaps. It provides two independent surfaces:

1. **Procedural 1D routines** — init, arithmetic, linspace, and cyclic shift on `ARRAY[*] OF REAL` or `ARRAY[*] OF STRING`. These compile into a single `matrix.pc` program and are imported via `%include matrix.klh`.

2. **2D matrix class template (`2dArray.klc`)** — a GPP-expanded class that generates a full-featured matrix program for any given dimensions and element type (REAL or INTEGER). Each instantiation becomes its own `.pc` file. The 4×4 variant is the foundation of `lib/pose/matpose`, which handles all homogeneous coordinate transforms in the 5-axis slicer.

Reach for this module when you need to:
- Initialize, transform, or sample arrays of reals
- Multiply, invert, or decompose matrices of fixed dimensions
- Build custom matrix types in your own module using the same class template

## Files

| File | Purpose |
|------|---------|
| `src/matrix.kl` | All `matrix__` 1D routines — implementation |
| `include/matrix.klh` | Public declarations + short aliases for 1D routines |
| `include/matrix.klt` | `OP_*` preprocessor constants (internal dispatch codes) |
| `include/matrix.private.klh` | Private declarations for `operations_1D` and `init_1D` |
| `include/classes/2dArray.klc` | GPP class template — expands into a named Karel program |
| `include/classes/2dArray.klh` | Class header — `declare_member` for all public methods |
| `include/classes/2dArray.private.klh` | Private method declarations for the class implementation |
| `include/classes/templates/carr3.klt` | Config: 3×3 REAL (rotation matrices) |
| `include/classes/templates/carr4.klt` | Config: 4×4 REAL (homogeneous transforms) |
| `include/classes/templates/carr10.klt` | Config: 10×10 REAL (general) |
| `include/classes/templates/carr23.klt` | Config: 2×3 INTEGER (general) |
| `include/classes/templates/carr305.klt` | Config: 30×5 REAL (path data) |
| `test/test_matrix.kl` | KUnit tests for all 1D routines |
| `test/test_matrix2.kl` | KUnit tests for all five built-in 2D class variants |

## API Reference

### 1D Routines

Import with `%include matrix.klh`. All routines live in `matrix.pc`.

#### Initialization

```
-- Fill arr with a constant
matrix__full_1D(arr : ARRAY[*] OF REAL; number : REAL)

-- Common fill shortcuts
matrix__zeros_1D(arr : ARRAY[*] OF REAL)   -- all 0.0
matrix__ones_1D(arr : ARRAY[*] OF REAL)    -- all 1.0
matrix__random_1D(arr : ARRAY[*] OF REAL)  -- random [0, 1], re-seeded per element

-- Evenly-spaced values between start and stop
-- endpoint=TRUE  -> [start, ..., stop] inclusive (N values over N-1 intervals)
-- endpoint=FALSE -> [start, ...) exclusive (N values over N intervals, NumPy default)
matrix__linspace(start_idx, stop_idx : REAL; arr : ARRAY[*] OF REAL; endpoint : BOOLEAN)

-- Write identity matrix into a native 2D array.
-- row is a dummy 1D array sized to the column count (Karel limitation workaround).
matrix__eye(arr : ARRAY[*,*] OF REAL; row : ARRAY[*] OF REAL)
```

#### Element-wise Arithmetic

All arithmetic routines require `arr1`, `arr2`, and `out_arr` to be the same length. UNINIT elements are passed through from whichever operand is initialized. Division by zero is silently skipped (`out_arr[i] = arr1[i]`).

```
matrix__add_1D(arr1, arr2 : ARRAY[*] OF REAL; out_arr : ARRAY[*] OF REAL)
matrix__sub_1D(arr1, arr2 : ARRAY[*] OF REAL; out_arr : ARRAY[*] OF REAL)
matrix__mul_1D(arr1, arr2 : ARRAY[*] OF REAL; out_arr : ARRAY[*] OF REAL)
matrix__div_1D(arr1, arr2 : ARRAY[*] OF REAL; out_arr : ARRAY[*] OF REAL)
```

#### Cyclic Shift

```
-- Cyclically rotate STRING or REAL array by `amount` positions.
-- updown=TRUE  -> shift up   (element at index i moves to i+amount, wrapping)
-- updown=FALSE -> shift down (element at index i moves to i-amount, wrapping)
-- unint_arr is a same-sized temp buffer -- declare it in the caller's VAR block.
-- The result is copied back into arr; unint_arr is modified as a side effect.
matrix__shift_s(arr : ARRAY[*] OF STRING; amount : INTEGER; updown : BOOLEAN; unint_arr : ARRAY[*] OF STRING)
matrix__shift_r(arr : ARRAY[*] OF REAL;   amount : INTEGER; updown : BOOLEAN; unint_arr : ARRAY[*] OF REAL)
```

#### Short Aliases

Every 1D routine has a 4-character `mtrx__` alias: `mtrx__full1R`, `mtrx__zros1R`, `mtrx__ones1R`, `mtrx__rand1R`, `mtrx__lspce_r`, `mtrx__eye`, `mtrx__shifts`, `mtrx__shiftr`, `mtrx__add1R`, `mtrx__sub1R`, `mtrx__mul1R`, `mtrx__div1R`.

---

### 2D Matrix Class

The class template is instantiated once per matrix shape you need. The result is a dedicated Karel program (e.g. `carr4.pc`) with methods named `carr4__zeros`, `carr4__mult`, etc.

**`.klt` config keys:**

| Key | Meaning | Example |
|-----|---------|---------|
| `type_name` | Name of the generated Karel type | `t_carr4` |
| `t_rows` | Number of rows | `4` |
| `t_columns` | Number of columns | `4` |
| `inner_type` | Element type -- `REAL` or `INTEGER` only | `REAL` |
| `UNIT_TESTING` | If defined, adds `test_row` and `testing` methods | (define or omit) |

**Type aliases generated by `define_type.m`:**
- `arry(type_name)` -- the class's 2D storage type (array of rows)
- `arrx(type_name)` -- the class's 1D row type

#### Creation

```
-- Parse a delimited string into a matrix.
-- row_delim splits rows; col_delim splits columns within each row.
carr3__create_array_from_string(sarr : STRING; col_delim : STRING; row_delim : STRING) : arry(type_name)

-- Parse a single row from a delimited string.
carr3__create_row_from_string(sarr : STRING; col_delim : STRING) : arrx(type_name)
```

#### Clearing and Setting

```
carr4__clear(arr : arry(type_name))
carr4__set_row(row_idx : INTEGER; row : arrx(type_name); out_arr : arry(type_name))
carr4__set_from_array(in_arr : ARRAY[*,*] OF REAL; out_arr : arry(type_name))
carr4__set_row_from_array(in_row : ARRAY[*] OF REAL; row_idx : INTEGER; out_arr : arry(type_name))
```

#### Getting

```
-- Square matrices only. Transposes internally to extract a column as a row.
carr4__get_col(col_idx : INTEGER; arr : arry(type_name)) : arrx(type_name)
```

#### Conversion

```
-- Copy class array to a native Karel ARRAY[*,*].
-- Only works reliably for square matrices (see Common Mistakes).
carr4__convert_to_array(in_arr : arry(type_name); out_arr : ARRAY[*,*] OF REAL)
```

#### Initialization

```
carr4__full(number : REAL; out_arr : arry(type_name))
carr4__zeros(out_arr : arry(type_name))
carr4__ones(out_arr : arry(type_name))
carr4__eye(out_arr : arry(type_name))
carr4__random(out_arr : arry(type_name))
```

#### Operations

```
-- Element-wise (available for all matrix shapes)
result = carr4__add(arr1, arr2)
result = carr4__subtract(arr1, arr2)
result = carr4__smult(arr, val)       -- scalar multiply: all elements * val

-- Square matrices only
result = carr4__mult(arr1, arr2)      -- matrix product
result = carr4__transpose(arr)
result = carr4__inverse(arr)          -- (1/det) * transpose(cofactor); aborts if |det| < 0.05
d      = carr4__det(arr)              -- scalar determinant (recursive cofactor expansion)
result = carr4__cofactor(arr)
t      = carr4__trace(arr)            -- sum of diagonal
```

#### Testing (requires `UNIT_TESTING` in `.klt`)

```
ok = carr4__test_row(expected, actual)   -- fuzzy compare one row (EPSILON = 0.005)
ok = carr4__testing(expected, actual)    -- compare entire matrix row by row
```

---

## Common Patterns

### 1. Using a built-in 2D template

```
-- In your .kl program:

%include define_type.m
%include carr4.klt
t_arr2d_ref(t_rows, t_columns, type_name, inner_type)

%class carr4('2dArray.klc','2dArray.klh','carr4.klt')

VAR
  mat : arry(type_name)   -- ARRAY[4] OF t_carr4
  row : arrx(type_name)   -- ARRAY[4] OF REAL

BEGIN
  -- X-rotation matrix
  carr4__eye(mat)
  row[1]=1 ; row[2]=0        ; row[3]=0         ; row[4]=0
  carr4__set_row(1, row, mat)
  row[1]=0 ; row[2]=COS(ang) ; row[3]=-SIN(ang) ; row[4]=0
  carr4__set_row(2, row, mat)
  row[1]=0 ; row[2]=SIN(ang) ; row[3]=COS(ang)  ; row[4]=0
  carr4__set_row(3, row, mat)
  row[1]=0 ; row[2]=0        ; row[3]=0         ; row[4]=1
  carr4__set_row(4, row, mat)
```

This is the exact pattern used by `lib/pose/matpose/matpose.kl` for `matpose__rotx/roty/rotz`.

### 2. Defining a custom matrix type in your module

When the built-in templates don't fit, create your own `.klt`:

```
-- include/mymods/mymat.klt
%include define_type.m
%defeval type_name t_mymat
%defeval t_rows 4
%defeval t_columns 4
%defeval inner_type REAL
-- omit %define UNIT_TESTING for production builds
```

Then in your `.kl`:

```
%include mymat.klt
t_arr2d_ref(t_rows, t_columns, type_name, inner_type)
%class cmat('2dArray.klc','2dArray.klh','mymat.klt')
```

This is exactly how `lib/pose/matpose` works -- it creates a `cmat` program from a custom `matarr.klt` (4×4 REAL, no UNIT_TESTING in production).

### 3. Parsing matrix data from a string

```
VAR
  mat : arry(type_name)
BEGIN
  -- semicolon-delimited rows, comma-delimited columns
  mat = carr3__create_array_from_string('2,-3,9;2,0,-1;1,4,5', ',', ';')
  -- mat[1] = [2,-3,9], mat[2] = [2,0,-1], mat[3] = [1,4,5]
```

### 4. linspace for trajectory sampling

```
VAR
  angles : ARRAY[10] OF REAL
BEGIN
  -- 10 angles from 0 to 360, inclusive
  matrix__linspace(0.0, 360.0, angles, TRUE)
  -- [0, 40, 80, 120, 160, 200, 240, 280, 320, 360]

  -- 5 values from 2.0 to 3.0, endpoint excluded (NumPy-style)
  matrix__linspace(2.0, 3.0, pts, FALSE)
  -- [2.0, 2.2, 2.4, 2.6, 2.8]
```

### 5. Cyclic shift for sliding-window buffers

```
VAR
  buf : ARRAY[5] OF REAL
  tmp : ARRAY[5] OF REAL   -- required temp buffer, same size as buf
BEGIN
  buf[1]=10 ; buf[2]=20 ; buf[3]=30 ; buf[4]=0 ; buf[5]=0
  -- shift up by 1: [10,20,30,0,0] -> [20,30,0,0,10]
  matrix__shift_r(buf, 1, TRUE, tmp)
  -- shift down by 3: [20,30,0,0,10] -> [0,10,20,30,0]
  matrix__shift_r(buf, 3, FALSE, tmp)
```

---

## Common Mistakes

| Mistake | Symptom | Fix |
|---------|---------|-----|
| Calling `mult`/`det`/`inverse`/`transpose`/`get_col`/`trace`/`cofactor` on a non-square class | Compile error: routine not defined | These are only compiled when `t_rows = t_columns`; use a square template |
| Omitting `unint_arr` temp buffer for `shift_s` / `shift_r` | Wrong result or runtime error | Declare a same-type, same-length array and pass it as the 4th argument |
| Using `matrix__eye` (1D API) with a class-typed variable | Karel type mismatch at compile | Use `carr4__eye(mat)` for class arrays; `matrix__eye` targets native `ARRAY[*,*]` |
| Setting `inner_type = STRING` in a custom `.klt` | GPP `%error` at compile time | Only `REAL` and `INTEGER` are valid inner types |
| Inverting a nearly-singular matrix | ABORT: "Matrix is not invertable" | Singularity threshold is `abs(det) < 0.05` (large); ensure matrix is well-conditioned |
| Forgetting `t_arr2d_ref(...)` after `%include yourmat.klt` | Type `t_carr*` or `t_mymat` undefined | `t_arr2d_ref` emits the type declaration; it must appear in the program body |
| `convert_to_array` on a non-square matrix | Spurious `ARR_LEN_MISMATCH` abort | Bug: both dimension checks use first-dimension length; only reliable for square matrices |

---

## Build Flow

`matrix` is Layer 2 and depends on `math`, `errors`, and `ktransw-macros`.

```
# From lib/matrix:
rossum .. -w -o        # resolve deps, generate build.ninja
ninja                  # compiles matrix.pc + carr3.pc, carr4.pc, carr10.pc, carr23.pc, carr305.pc
kpush                  # deploy to controller

# With tests:
rossum .. -w -o -t     # adds test_matrix.pc, test_matrix2.pc (pulls in KUnit, Strings)
ninja && kpush
# http://<robot-ip>/KAREL/kunit?filenames=test_matrix,test_matrix2
```

Each `%class` expansion compiles to its own `.pc`. When another module creates a custom class (e.g. `cmat`), that module's build will produce `cmat.pc` separately.

See the [Ka-Boost readme](../../readme.md) for full toolchain setup and build instructions.
