# Feature spec: loop_link, a Python link to geoviewer

**Feature**: `006-loop-link` | **Created**: 2026-11-05 | **Status**: Draft

**Input**: A Python library, `packages/loop_link`, that sends 3D objects with
their properties to geoviewer, receives them back, and shares the selection
(the highlight) between Python and the viewer. The LoopStructural QGIS plugin
uses the same library.

<!--
Relation to other specs:
- Depends on geoviewer spec 066-selection-and-fetch (selection messages and a
  way to fetch objects from the viewer). That spec is in the geoviewer
  repository. A tool selection opens the property table of one object (nodes
  or cells), so a selection holds one object. A fetch of an object over
  256 MiB arrives in parts. Stories 2, 3 and 4 wait for it, and so does the
  map-name list of story 7 (066, FR-009). Stories 1, 5, 6 and the rest of
  story 7 work with the viewer as it is.
- The viewer lets one data client send objects at a time (geoviewer spec 003).
  Tool connections (up to 8) can select, read, fetch, filter and colour.
- Replaces the push-only client in the plugin repository
  (`loopstructural/main/geoviewer`). Decision of 2026-10-31: the plugin uses
  this library directly.
-->

## User scenarios and tests

### User story 1: Send objects to the viewer (priority: P1)

A Python user sends points, vectors, lines, surfaces, structured grids, tetra
meshes and drillholes to a running viewer. Each object keeps its properties.
The user sees the objects in the viewer.

**Why this priority**: All other stories need objects in the viewer.

**Independent test**: Start the fake viewer (FR-015). Send one object of each
type with one property. Read the inventory of the fake viewer and check the
ids and types. A manual check with the real headless viewer repeats the test
(CI cannot run it).

**Acceptance scenarios**:

1. **Given** a running viewer and a saved token, **when** the user sends an
   object, **then** the library finds the viewer with no port or token entry
   and the object is in the viewer scene.
2. **Given** an object that was sent, **when** the user sends the same id
   again, **then** the viewer replaces the object and the revision increases.
3. **Given** no running viewer, **when** the user sends an object, **then** the
   library gives a clear error and does not start a viewer.
4. **Given** an object with properties on its points or cells, **when** it is
   sent, **then** the viewer shows the same property names and values.
5. **Given** one object of each of the seven types, **when** the user sends
   them, **then** the viewer holds seven objects with the correct types.
6. **Given** an x, y, z direction vector, **when** it is sent, **then** the
   viewer holds the equal orientation or lineation.
7. **Given** one Python program that holds the data connection of a viewer,
   **when** a second program tries to send, **then** the second program gets
   the library error `busy` and the objects of the first program do not
   change. The second program can still select, read, fetch, filter and
   colour on the tool connection.
   **Given** the first program releases the data connection, **when** the
   second program sends, **then** the viewer accepts it and keeps the objects
   of the first program.
8. **Given** an install of the core only, **when** the library is imported,
   **then** it needs neither Qt nor the optional extras.

### User story 2: Select a cell from Python (priority: P1)

The user sends an object id and one or more cell indices. The viewer selects
them. The selection is the same selection that a user makes in the viewer.

**Why this priority**: It is the main way to link a result in Python to the
picture.

**Independent test**: Send a structured grid. Select cell 12 from Python.
Read the selection of the viewer and check that it holds cell 12.

**Acceptance scenarios**:

1. **Given** an object in the viewer, **when** Python selects cell indices,
   **then** the viewer outlines those cells.
2. **Given** an index that is outside the object, **when** Python selects it,
   **then** the library gets an error and the selection does not change.
3. **Given** an object that was replaced after Python read its indices,
   **when** Python selects with the old revision, **then** the viewer refuses
   the selection.
4. **Given** a viewer without the selection messages, **when** Python selects,
   **then** the library gives an error that names the missing feature.

### User story 3: Read the selection in Python (priority: P1)

The user selects objects or cells in the viewer. Python reads the selection on
request. The library can also poll the selection and call a function when it
changes. The viewer does not push a selection change.

**Why this priority**: It closes the loop from the picture to Python.

**Independent test**: A second client selects a cell in the viewer. Python
reads the selection and gets the object id and the cell index.

**Acceptance scenarios**:

1. **Given** a selection in the viewer, **when** Python asks for it, **then**
   Python gets each selected object id, the kind of location (cell, vertex)
   and the indices.
2. **Given** a polling callback in Python with a set interval, **when** the user
   changes the selection in the viewer, **then** Python calls the callback
   once for the change, at the next poll.
3. **Given** a selection of more than 100,000 cells, **when** Python asks for
   it, **then** the answer arrives as an array and not as text.

### User story 4: Receive objects from the viewer (priority: P2)

Python gets an object from the viewer as the same kind of object that it can
send: geometry and properties. The objects include those that the viewer made
by an import, for example a drillhole import or a file.

**Why this priority**: It makes the link two-way. It needs a new viewer
feature, so it comes after send and select.

**Independent test**: Send an object, fetch it back, and compare geometry and
properties with the original.

**Acceptance scenarios**:

1. **Given** an object in the viewer, **when** Python fetches it, **then**
   Python gets the geometry and all properties with the same values.
2. **Given** an object type, **when** Python sends it and fetches it back,
   **then** the result equals the original (a round trip).
3. **Given** a request for an id that does not exist, **when** Python fetches
   it, **then** the library gives a clear error.

### User story 5: The QGIS plugin uses the library (priority: P2)

The plugin sends model objects and reads the selection through this library.
The library does not need Qt. The plugin runs it on a worker thread.

**Why this priority**: It removes the second copy of the client code.

**Independent test**: Push the same model with the library and with the old
client to the fake viewer, and compare the messages that the fake viewer
receives.

**Acceptance scenarios**:

1. **Given** the plugin, **when** it pushes a model, **then** it calls only the
   library and the objects in the viewer are the same as with the old client.
2. **Given** a slow viewer, **when** the plugin sends, **then** the QGIS user
   interface does not block.

### User story 6: Filter objects from Python (priority: P2)

The user sets a data filter on an object in the viewer from Python. The filter
is the filter of the viewer (spec 023 of geoviewer): a range filter on a number
property, or a value filter on a category property. The viewer shows only the
elements that pass.

**Why this priority**: The user finds a part of a model, for example the cells
with a grade above 0.5, without a manual step in the viewer.

**Independent test**: Send a grid with a number property. Set a range filter.
Check that the viewer reports the filter as on, with the count of elements
that pass. Use the fake viewer in CI.

**Acceptance scenarios**:

1. **Given** an object with a number property, **when** the user sets a range
   (minimum, maximum, and keep or drop NaN), **then** the viewer shows only the
   elements in the range.
2. **Given** an object with a category property, **when** the user sets the
   values to show, with or without the exclude switch, **then** the viewer
   shows only the matching elements.
3. **Given** a filter, **when** the user turns it off or removes it, **then**
   the viewer shows all elements again. The data does not change.
4. **Given** a property that does not exist, or a property of a type that has
   no filter (text, or one value for the whole object), **when** the user sets a
   filter, **then** the library gives a clear error.
5. **Given** a filter that takes time to build, **when** the user sets it,
   **then** the call can wait for the end, or return at once, and the user can
   cancel the build.
6. **Given** a new revision of an object, **when** the user sends it, **then**
   the filters of the object stay and apply to the new data.

### User story 7: Set the colouring of an object from Python (priority: P2)

The user chooses which property colours an object, and sets the colour map,
the range and the scale. The user can also set the plain colour of an object.
These are the colour settings of the viewer (geoviewer specs 018 and 047).

**Why this priority**: A script can show a result at once with a useful colour,
for example a grade with a log scale, with no manual step in the viewer.

**Independent test**: Send a grid with a number property. Set the property as
the colouring, the map `coolwarm`, a range and the log scale. Check the
settings with the fake viewer.

**Acceptance scenarios**:

1. **Given** an object with properties, **when** the user chooses a property
   to colour by, **then** the viewer colours the object by that property. The
   user can also choose no property.
2. **Given** a number property, **when** the user sets a colour map by name,
   **then** the viewer uses that map. The names are those of the viewer.
3. **Given** a number property, **when** the user sets a range, or clears it,
   **then** the viewer uses that range, or the range of the data.
4. **Given** a number property, **when** the user chooses the linear or the log
   scale, **then** the viewer uses that scale.
5. **Given** a property that is a category or a boolean, **when** the user sets
   a map, a range or a scale, **then** the library gives a clear error,
   because these properties have none.
6. **Given** a map name that the viewer does not know, **when** the user sets
   it, **then** the library gives a clear error that lists the names.
7. **Given** the user sends an object with the colour settings, **when** the
   viewer receives it, **then** the settings apply to the object. The
   settings go in the send, or in a follow-up call that the library sends at
   once [NEEDS CLARIFICATION: the plan checks if the send message can carry
   the settings; the user-visible result is the same].

### Edge cases

- The token is missing or wrong. The library says how to set the token.
- The viewer is an older version without the selection messages. Sending works.
  Select and read give a clear error that names the missing feature.
- An object is larger than one frame (256 MiB) or larger than the object limit
  (2 GiB). The library refuses it with an error before it sends.
- The origin of the viewer is different from the origin of the client.
- The viewer closes in the middle of a send.
- The user has a property that the protocol cannot hold (for example an
  object array).
- Two Python programs connect to the same viewer. Only one holds the data
  connection (story 1, scenario 7).
- A filter or a colour command changes only how the viewer shows an object.
  It is not an edit, and it has no owner check or undo step.

## Requirements

### Functional requirements

- **FR-001**: The library must find a running viewer and its token with no
  manual entry.
- **FR-002**: The library must send these object types with their properties:
  points, vectors, lines, surfaces, structured grids, tetra meshes and
  drillholes.
- **FR-003**: The library must keep the revision of each object it sent. It
  must replace an object that has the same id, and the revision must
  increase.
- **FR-004**: The library must select objects and cells in the viewer by
  object id and index.
- **FR-005**: The library must read the selection of the viewer on request,
  and must offer a poll that calls a function when the selection changes. The
  user sets the poll interval (default 0.5 s). The poll runs only when the
  user starts it.
- **FR-006**: The library must fetch any object of the scene from the viewer,
  also an object that another client or an import made, and return the
  geometry and the properties. A fetch does not change the scene. Coordinates
  are returned as world coordinates. Each array keeps the data type that the
  viewer holds.
- **FR-014**: The library must hold at most one data connection for each
  program. When the viewer refuses a second data client with `busy`, the
  library must give the error `busy` and must not change the objects of the
  first program. A program without the data connection can still select, read
  the selection, fetch, filter and colour. The library must let a program
  release the data connection and keep the tool connection, so that another
  program can send. The viewer keeps the objects that the first program sent.
- **FR-007**: The library must not need Qt. A plugin adapter handles threads.
- **FR-008**: The library must give a library error for each failure in the
  edge cases. Each failure has its own error type, with a text that says what
  to do.
- **FR-009**: The library must keep an API that the plugin can use in place of
  its own client. The plugin change, and the removal of the old client, are
  done in the plugin repository after the library ships.
- **FR-010**: The library must accept these inputs: pyvista meshes and grids
  and plain numpy arrays with a dict of properties, for all seven object
  types; LoopStructural objects for structured grids and surfaces; pandas
  tables (collar, survey, logs) for drillholes only. Any other combination
  gives a clear error that names the accepted inputs.
- **FR-011**: The library must accept a vector as an orientation (dip
  direction and dip), a lineation (trend and plunge) or an x, y, z direction.
  It converts an x, y, z direction before it sends.

- **FR-012**: The library must own the only Python codec of the protocol. The
  core needs numpy and a WebSocket client only. The adapters for pyvista,
  LoopStructural and pandas are optional extras.
- **FR-016**: The change must add the package to `ROADMAP.md` and to the
  per-package release configuration.
- **FR-017**: The package must support Python 3.10 and later. This is the
  minimum of the package. It is higher than the lint target of the repo
  (py39), so the plan must set the target for this package.
- **FR-018**: The library must set, turn on or off, and remove a range filter
  and a value filter (with the exclude switch) on an object in the viewer, by
  the name of a property. The library uses the existing filter command of the
  viewer. This needs no change in the viewer.
- **FR-019**: The library must let the user wait for a filter build to end, or
  continue at once, and must cancel a filter build on request.
- **FR-020**: The library must set the property that colours an object, the
  colour map by name, the range (and clear it), and the linear or log scale.
  It must also set the plain colour of an object. The library uses the
  existing colour commands of the viewer. It must get the list of map names from
  the viewer (the command list of the viewer), not from a fixed list in the
  library. This needs a small change in the viewer (geoviewer spec 066, FR-009).
- **FR-021**: The library must let the user give the colour settings in the
  call that sends an object. The library applies them with the send or with a
  follow-up call, so that the user needs one call.
- **FR-015**: The package must include a fake viewer: a test server that
  speaks the protocol. CI uses it, because CI cannot run the real viewer.
- **FR-013**: Stories 2 to 4 need geoviewer spec 066 in the viewer. The
  library must detect a viewer without it, and must give a library error that
  names the missing feature for these calls.

### Protocol and data

- **Messages**: uses the messages of protocol 1.14 for send, and the existing
  commands for filter and colour. Uses the selection messages and the fetch
  message of geoviewer spec 066: set the selection, read the selection, and
  fetch an object. The exact names and the protocol version come from that
  spec. The viewer sends no selection change message.
- **Objects and data**: object types and data types of the viewer. A "vector"
  is a Points object with an orientation or lineation property group. An
  x, y, z direction is converted to dip direction and dip, or to trend and
  plunge, before it is sent. A fetch returns the angles that the viewer
  stores, and a helper converts them to x, y, z (assumption, 2026-11-05: the
  viewer keeps the angles, so the round trip of an x, y, z input is equal
  only up to the conversion).

### Key entities

- **Object**: one item in the viewer scene, with an id, a type, geometry and
  properties.
- **Property**: a named array of values on the points, cells or the whole
  object.
- **Viewer connection**: the link from Python to one viewer. It has a tool
  connection and, when the program sends objects, the one data connection. It
  holds the revision of each object it sent.
- **Library error**: an error type of the library with a text. The types
  include no viewer, wrong token, `busy`, missing feature, invalid input and
  object too large.
- **Selection**: a list of object ids, and for each object a location kind and
  an array of indices.

## Success criteria

- **SC-001**: A grid of 1,000,000 cells with one float property is sent and
  appears in the viewer inventory in less than 5 s on the same computer.
- **SC-002**: A send followed by a fetch gives geometry and properties that
  are equal for all seven object types, within a relative tolerance of 1e-6
  (the precision of `f32`).
- **SC-003**: A selection of 100,000 cells goes from Python to the viewer in
  one message.
- **SC-004**: The fake viewer reads all 100% of the frames in the geoviewer
  folder `fixtures/protocol/python`, and the library writes frames that the
  geoviewer tests accept.
- **SC-005**: The package installs and imports on Python 3.10 and later, with
  0 imports of Qt in the core.
- **SC-006**: `ROADMAP.md` and the release configuration list `loop_link`
  (2 entries).
- **SC-007**: A colour command, or a filter command, gets its reply from the
  fake viewer in less than 250 ms for an object with 1,000,000 cells. The time
  of a filter build is not counted.

## Out of scope

- The geoviewer changes (geoviewer spec 066). They are a separate spec.
- Several programs that send objects at the same time. The viewer allows one
  data client.
- A separate highlight that is not the user selection.
- Starting or installing the viewer.
- Removing the old client from the plugin repository.
- A push message for selection changes.
- A custom colour map (own colours or stops). Only the named maps of the
  viewer.
- Reading the colour settings back from the viewer.
- Reading the filters back from the viewer. The viewer has no command for it.
- Filters on drillhole log sets (the filter of geoviewer spec 049).
- Editing an object that another client owns.
- A link to the browser (wasm) build of the viewer. The browser cannot run
  the server that the library connects to. A later spec can add a notebook
  widget or a Python-side server.

## Clarifications

- 2026-11-05: Q: Where does the library live? A: `packages/loop_link`.
- 2026-11-05: Q: Highlight and selection: same or separate? A: Same. The
  highlight is the user selection.
- 2026-11-05: Q: Which object types? A: Points, vectors, lines, surfaces,
  structured grids, tetra meshes and drillholes, all with properties. Send and
  receive.
- 2026-11-05: Q: Name of the spec? A: `006-loop-link`.
- 2026-11-05: Q: What is a vector? A: Both: orientation or lineation, and an
  x, y, z direction (converted before sending).
- 2026-11-05: Q: Which Python inputs? A: pyvista meshes and grids,
  LoopStructural objects, plain numpy arrays, pandas drillhole tables.
- 2026-11-05: Q: Who owns the Python protocol codec? A: `loop_link` alone. It
  is tested against the geoviewer frames in `fixtures/protocol/python`. The
  plugin drops its copy.
- 2026-11-05: Q: How are the adapters packaged? A: Optional extras (pyvista,
  LoopStructural, pandas). The core needs numpy and a WebSocket client only.
- 2026-11-05: Q: Viewer side of stories 2 to 4? A: A new geoviewer spec for
  protocol 1.15. Stories 2 to 4 wait for it. Stories 1 and 5 work with 1.14.
- 2026-11-05: Q: Speed of a send? A: Under 5 s for 1,000,000 cells.
- 2026-11-05: Q: Which objects can be fetched? A: Any object of the scene.
- 2026-11-05: Q: Selection events? A: Polling only. The viewer pushes no
  selection message; the library polls. This lowers the viewer work to
  `selection.set` and `selection.get`.
- 2026-11-05: Q: Many clients? A: Each client owns its ids. The selection is
  shared. (Replaced in the second review: the viewer allows one data client.)
- 2026-11-05 (second review): Q: Can two programs send at once? A: No. The
  viewer gives `busy` to a second data client (geoviewer spec 003). The second
  program can use the tool connection.
- 2026-11-06: Q: Should the viewer allow several data clients? A: No. Keep one
  data client. The library lets a program release the data connection so that
  another program can send (FR-014).
- 2026-11-05 (review): Q: Tolerance of a round trip? A: Relative 1e-6.
- 2026-11-05 (review): Q: Data type of a fetch? A: The type the viewer holds;
  coordinates are world coordinates.
- 2026-11-05 (review): Q: Which inputs for which types? A: pyvista and numpy
  for all seven; LoopStructural for grids and surfaces; pandas for drillholes.
- 2026-11-05 (review): Q: Poll interval? A: Set by the user, default 0.5 s.
- 2026-11-05 (review): Q: Minimum Python? A: 3.10 and later (FR-017).
- 2026-11-05 (review): Q: ROADMAP and release entries? A: Required in this
  spec (FR-016).
- 2026-11-05 (review): Q: Geoviewer spec not yet written? A: Keep the
  dependency as text. The plan covers stories 1 and 5 first. Fill in the
  number when the geoviewer spec exists.
- 2026-11-05: Q: Can Python also send filters? A: Yes. The viewer already has
  the command `filter.edit` on the tool connection (range, values, exclude,
  on or off, remove, cancel). Story 6, FR-018 and FR-019. No viewer change.
- 2026-11-05: Q: Can Python change the colour map? A: Yes. The viewer has the
  commands `colour.set_active`, `colour.set_map`, `colour.set_range`,
  `colour.clear_range`, `colour.set_scale` and `object.set_colour`. See story 7
  and the requirements for colour. The viewer has three named maps (`viridis`, `grey`,
  `coolwarm`). A custom map is not in scope.
- 2026-11-05: Q: Where do the colour map names come from? A: From the viewer.
  Today the input `name` of `colour.set_map` is free text and the command list
  does not give the names, so spec 066 adds the names to the schema.
- 2026-11-05: Q: Old plugin client? A: Removed in the plugin repository after
  the library ships, not in this spec.
