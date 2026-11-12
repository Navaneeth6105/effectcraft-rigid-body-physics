# 2D Physics Panel for EffectCraft

A script panel that makes shape and text layers behave as rigid bodies. It reads each layer's real outline, simulates gravity, collisions, bounce, friction and spin, then bakes the result into ordinary Position and Rotation keyframes (editable, fast to render, no expressions).

## Install and run

1. Put PhysicsRig.jsx in EffectCraft's scripts folder, or run it from the Scripts menu (File \> Scripts \> Run Script File), or paste it into the Script Console.
2. A "2D Physics" panel appears (docked if the app opens it as a panel, otherwise a floating palette).

## Basic use

1. Select the layers to simulate and put the playhead where the simulation should start.
2. Tick **Gravity enabled**, adjust the values, press **Simulate selected layers**.
3. Change something and press it again to re-bake. Keys from the playhead onward are replaced.
4. **Clear keys from current time** removes baked keys. Undo also works.

| Setting | What it does |
| --- | --- |
| Gravity / Gravity starts at | Strength in px/s², and the time (seconds after the start) when it switches on |
| Default bounce, Default friction | 0 to 1. Used by every layer that has no Bounce/Friction controls of its own |
| Initial push | Random starting velocity and spin |
| Duration | Seconds to bake |
| Walls | Floor, left, right, top: the composition edges |
| Collider | **Follow shape outline (concave)**, **Convex hull**, or **Bounding box** |
| Carry velocity from keyframes | A layer you animated keeps its speed and spin at the playhead when physics takes over |

## Layer names decide what a layer is

The prefix must be followed by a **space, underscore or dot**. Anything else is just a normal layer.

| Name | Behaviour |
| --- | --- |
| Static\_x, Static x, Static.x, and the same with Fixed | Immovable obstacle (any number of them; they need not be selected) |
| Anim\_x, Anim x, Anim.x, and the same with Keyed or Kinematic | Obstacle that follows its own position/rotation keyframes and pushes other objects |
| any other name, selected | Moved by physics |

So Animal, Staticfloor or Fixedpost are normal layers, and Anim\_arm or Static floor are obstacles. Not case sensitive.

## Per-object bounce, friction and density

Press **Add Bounce/Friction/Density controls** with layers selected. Each gets three Slider Controls named Bounce, Friction and Density (filled with the current panel defaults). Change them per layer; select several layers to edit them together. The script reads them at the playhead when you press Simulate.

- Bounce: 0 = dead, 1 = very bouncy. Friction: 0 = ice, 1 = grippy. Density: heavier layers push lighter ones.
- Two touching bodies combine as: bounce = the higher of the two, friction = the geometric mean.
- Floor and walls adopt the object's own values. A Static/Anim layer can have its own controls too (a rubber floor, an icy ramp).
- If your build can't add effects by script, add Slider Controls by hand and rename them exactly Bounce, Friction, Density.

## Letters and separate shapes

Text is one layer, so it moves as one block. Press **Split selected into letters / shapes**:

- **Text layer** becomes one text layer per letter (spaces skipped; multi-line text and left/centre/right alignment supported), placed exactly where the letters were. The original is hidden, the new letters are selected and given the material controls. Then press Simulate. Each letter collides as its bounding box.
- **Shape layer with several shapes** becomes one layer per top-level shape group. This is how you get true letter outlines: convert the text to shapes first (if EffectCraft has Layer \> Create Shapes from Text), then Split. Now "C", "U" and "L" collide by their real concave outline and can hook into each other.

## Obstacles that move, and keyframed hand-over

- Make a layer Anim\_paddle, keyframe it as usual, select only the layers that should fall, press Simulate.
- Animate a layer normally, put the playhead where it should be let go, select it, tick **Carry velocity**, Simulate.
- Opacity, colour, effects and text are never touched. Scale is read once at the start.

## How shapes are treated

- Shape layers: rectangles, ellipses, stars and Bezier paths (also inside groups) are read from the real outline. Concave shapes are split into convex pieces that move as one rigid body, so objects can sit inside a cup or bowl. "Convex hull" shrink-wraps instead.
- Text, footage, solids: the layer's bounds box.
- Rotation happens around the centre of mass, whatever the anchor point is.

## Limits

- Text letters use the bounding box. Kerning in the split is approximate.
- Shapes with holes are solid. Self-intersecting paths fall back to the convex hull.
- Detailed outlines are simplified to about 40 points. Heavy curves create many pieces and are slower.
- Animated obstacles don't collide with walls or Static layers; don't animate one so it squeezes an object into a wall.
- Unparented 2D layers only. Skipped layers are listed in the panel.
- Baked at the composition frame rate; physics stepped at 120 Hz internally.

## Testing status

The physics (bounce ratios, friction on a slope, density, stacking, cups, animated obstacles) and all layer handling (name rules, material sliders, splitting text and shape groups, baking) were tested against a mock of the After Effects scripting model, not inside EffectCraft itself. Splitting text relies on duplicating a layer, rewriting its text and measuring it; if EffectCraft doesn't support one of those, the panel says so and suggests converting to shapes. Send me any error text the panel shows.
