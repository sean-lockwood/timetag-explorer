# Selections in the `_tagsel.fits` file

**Status: implemented in package 1.7.0 (selection script 3.0); `CURVE` column added in 1.8.0 (selection script 3.1); planet-map columns added in 1.9.0 (selection script 3.2).** It replaces the single time window and single region
of selection script 2.x with any number of named time ranges and regions (selection script 3.0). Nobody depends on the
2.x layout, so there is no compatibility mode.

## In short

* Each time range or region owns one **bit** in a 16-bit flag column of the `EVENTS` table: `SEL_TIME` or `SEL_REGION`.
  An event inside several regions has several bits set; 0 means in none.
* Two table extensions, `TIMESEL` and `REGSEL`, define the bits: shape, frame, label, color, visibility.
* A region built from others (union, intersection, difference) is an **expression** over them and takes no bit.
* The existing `EDGE_OK` flag, positions, weights and processing keywords are unchanged.

## Flag columns

| Column | Type | Meaning |
|---|---|---|
| `SEL_TIME` | `I` with `TZERO = 32768` (unsigned 16-bit) | bit *n* set: the event is in time range with `BIT = n` |
| `SEL_REGION` | the same | bit *n* set: the event is in region with `BIT = n` |

Bit 0 is the least significant. astropy returns these as `uint16`; some other readers return signed values, in which
case bit 15 looks like a sign. At most 16 time ranges and 16 regions have a bit; derived rows do not count.
A column is written only if at least one selection of its kind exists.

## `TIMESEL` table (`EXTNAME = 'TIMESEL'`, `EXTVER = 1`)

| Column | Type | Meaning |
|---|---|---|
| `NAME` | string | short identifier (`T1`, `orbit2`): starts with a letter, then letters, digits, `_`; unique ignoring case; used in expressions |
| `NOTE` | string | the user's label, free text |
| `BIT`, `MASK` | `I`, `J` | bit number and its value `2**BIT`; `-1` and `0` for a derived row |
| `TSTART`, `TSTOP` | `D`, seconds | time range `TSTART <= TIME < TSTOP`, seconds since the exposure start (the `TIME` column); `NaN` for a derived row |
| `EXPR` | string | derived rows only: expression over earlier rows; empty otherwise |
| `COLOR` | string | `#rrggbb`, the color in plots |
| `VISIBLE` | `L` | shown in plots |

## `REGSEL` table (`EXTNAME = 'REGSEL'`, `EXTVER = 1`)

Same `NAME`, `NOTE`, `BIT`, `MASK`, `EXPR`, `COLOR`, `VISIBLE` columns, plus:

| Column | Type | Meaning |
|---|---|---|
| `FRAME` | string | reference frame the region was drawn in: `image`, `rest`, `rest_north`, `planet`, `planet_north` |
| `SHAPE` | string | `BOX`, `CIRCLE` or `POLYGON`; empty for a derived row |
| `PARAMS` | `nD`, pixel | shape numbers, in float64, padded with `NaN`; `n` is the longest row |
| `DESCR` | string | the one-line description the logs already use, for people |
| `CURVE` | `L` | the Explorer also draws this region's light curve. With `VISIBLE` (its overlay) it gives four display states: off, curve only, overlay only, on. Added in selection script 3.1; an older file has none and it reads as true. Nothing outside the Explorer needs it. |

`PARAMS` by shape:

| Shape | Numbers | Test |
|---|---|---|
| `BOX` | `x0, x1, y0, y1` | `x0 <= x < x1` and `y0 <= y < y1` |
| `CIRCLE` | `cx, cy, r` | `(x-cx)^2 + (y-cy)^2 <= r^2` |
| `POLYGON` | `x1, y1, x2, y2, ...` (closing edge implied) | inside by the even-odd rule; points exactly on an edge are not guaranteed either way |

A derived row has no `FRAME`; its frames are those of the rows in its expression. A new shape later (ellipse, annulus)
is a new `SHAPE` name with documented numbers; readers that do not know it can still use the bit column.

## Coordinates

* Regions use the **same 1-indexed hi-res pixel coordinates as the position columns** (`XCORR`, `XREST`, ...), so a
  region from the table can be applied by hand to those columns. The shift to 0-indexed arrays happens in one place in
  the package (the region is drawn at +1, plotted at −1, and one helper hands array-based code a 0-indexed copy).
* The flag is computed at each event's own position, on exactly the float32 values written to the position columns, so
  re-applying a region from the table by hand to those columns gives the same flags.
* Planet-map frames (`map_*`, **experimental**) are in **degrees** and are never shifted by one. A region drawn on the rectangular map is tested on the `MAPLON` (west longitude, unwrapped around `MAPLCEN`) and `MAPLAT` (planetographic latitude) columns; one drawn on a polar map, on `MAPNX`/`MAPNY` (north) or `MAPSX`/`MAPSY` (south), degrees of colatitude in the azimuthal-equidistant plane. Events that miss the planet are `NaN` there and are in no region. `MAPLON`/`MAPLAT` are written whenever the exposure is of Jupiter or Saturn (or one of their moons), the polar columns only when a region was drawn on them. The columns depend on the planet registration, which is recorded in the primary header (`MAPDX`, `MAPDY`: planet-centre offsets in hi-res pixels; `MAPALT`: emission altitude in km; `MAPLCEN`: the centre longitude of the rectangular map) and in the script's `CORRECTIONS`. Changing the registration changes these coordinates, so a region drawn under one registration is not valid under another.

## Leaving things out

In the package (not in the file, which always has every column) a selection needs only its numbers. A missing `NAME`
becomes `T1`, `T2`, ... or `R1`, `R2`, ... in table order (numbers already taken by a given name are skipped); a missing
bit becomes the lowest free one, around any bits given; a missing color comes from a default list; `VISIBLE` is true.
A derived row must have a name, since other rows refer to it. Give names or bits when the order may change.

## Choosing a selection

A time range or region is picked by its **`NAME`** (case-insensitive: `'Oval'` finds `oval`) or by its **bit position**
(`0`, `1`, `2`, ...; the bit number, not the mask value `1`, `2`, `4`, ...). A name never looks like a number, so the two
cannot be confused. A derived row has no bit and can only be picked by name.

## Expressions

`EXPR` combines `NAME`s (case-insensitive) of rows **above it** in the same table, so there are no cycles:

| Operator | Meaning |
|---|---|
| `A \| B` | union |
| `A & B` | intersection |
| `A ^ B` | in A or B but not both (exclusive or) |
| `A - B` | A and not B (same as `A & ~B`) |
| `~A` | not A |
| `( )` | grouping |

Precedence is `~`, then `&` and `-`, then `^`, then `\|` (as in Python); use parentheses to be clear. The text is parsed by the package, never
passed to `eval`. Time ranges and regions are separate tables and cannot refer to each other.
Combining across frames is exact, because each bit is already an event-level flag.

## What a script writes

* Position columns for the image frame (when corrections are on, or a region uses it) and for **every frame a region
  uses**, north-up or not, so each region can be re-applied by hand. `INCLUDE_NORTH_UP` still adds the rest.
* If any region's frame cannot be built (no SPT file, no Horizons answer), the script **stops before writing anything**.
* Header: `SELTMIN`, `SELTMAX`, `SELREG` and `SELFRAME` go away (the tables hold this); `TTETSEL` and `TTERSEL` stay,
  with the other `TTE*`, `RANDSEED`, `VIG*` and `SELEDGE*` keywords.
* `TIMESEL` and `REGSEL` follow the existing extensions, which are carried over untouched.

## Using it

The package reads the tables and hides the bit arithmetic:

```python
sel = read_selections('od8k57itq_tagsel.fits')            # time ranges and regions, with their numbers and colors
t = Table.read('od8k57itq_tagsel.fits', hdu=1)
keep = sel.in_time(t, 'orbit2') & sel.in_region(t, 'Oval') & t['EDGE_OK']   # by NAME, any case; derived names work too
keep = sel.in_time(t, 1) & sel.in_region(t, 0)             # or by bit position (the same selections here)
keep = (t['SEL_REGION'] & 0b101) != 0                      # or by hand: in the region with bit 0 or bit 2
```

Worked out when asked, not stored: `sel.time_exposure('orbit2', gti)` (good time inside the range),
`analysis.time_exposure(ds, times, name)`, and `analysis.region_stats(ds, regions, name, corr, times, time_name)` (area
on the Explorer's 4 x 4 pixel bins, and the mean, minimum and maximum exposure over them). A derived region must be made
of regions in one frame for its area and exposure.

Example tables for one exposure. `TIMESEL`:

| NAME | NOTE | BIT | TSTART | TSTOP | EXPR |
|---|---|---|---|---|---|
| `orbit1` | first orbit | 0 | 0 | 1500 | |
| `orbit2` | second orbit | 1 | 3000 | 4500 | |
| `both` | both orbits | −1 | NaN | NaN | `orbit1 \| orbit2` |

`REGSEL`:

| NAME | NOTE | BIT | FRAME | SHAPE | EXPR |
|---|---|---|---|---|---|
| `oval` | main oval | 0 | `planet_north` | `POLYGON` | |
| `io` | Io footprint | 1 | `rest` | `CIRCLE` | |
| `outside` | oval without footprint | −1 | | | `oval - io` |

(Other columns omitted.) "Brightness of the oval during orbit 2" is `in_time(t, 'orbit2') & in_region(t, 'oval')`. The file does not store such
pairings; the analysis picks them. (A `SELECTIONS` table that names them could be added later without changing any of
the above, for example to reload app state.)

## Decided

* `NAME` is the identifier; `NOTE` is the free label (ASCII; other characters become `?`). Selections are picked by name or bit position (above).
* 16 bits per column; a 32-bit word can follow later if needed.
* Precedence as in Expressions.
* Area and exposure are worked out when asked, not stored (exposure varies across a region and depends on the time range).
* A region's frame that cannot be built stops the script before anything is written; every frame a region uses gets its columns.
* Color and visibility are stored, for reloading a session later.
* **Map frames are left out for now**: a region drawn in a `map_*` frame makes the script stop. They will be added once
  the user has to determine and state an accurate planet registration.

## Not yet

* Reloading the Explorer's state from a file (the columns for it are there).
* Map-frame regions (above).
* More than 16 bits per column.
