# Cmgui visualization file, generated by opendihu at 16/5/2018 20:33:12
# run by cmgui laplace3d.com
$group = "Region"
$n = 1;

# read files
$fname = "laplace3d"
$nodefile = $fname . ".exnode"
$elemfile = $fname . ".exelem"
gfx read nodes node_offset $n $nodefile;
gfx read elements node_offset $n line_offset $n face_offset $n element_offset $n $elemfile;
$n+=1500;

# create materials
gfx create material muscle_transparent normal_mode ambient 0.4 0.14 0.11 diffuse 0.5 0.12 0.1 emission 0 0 0 specular 0.3 0.5 0.5 alpha 0.25 shininess 0.2
gfx create spectrum new log exaggeration 10 left range 0 1 extend_above extend_below rainbow colour_range 0 1 component 1;


# create faces
gfx define face egroup $group

# add line representation
gfx modify g_element "/" lines

# add surface representation with muscle material
gfx modify g_element $group surfaces material muscle_transparent
gfx modify g_element $group surfaces material muscle_transparent data solution spectrum new

# add streamlines
gfx modify g_element $group streamlines cell_centres coordinate geometry data gradient length 500000 circle_extrusion vector gradient travel_scalar

# add axes
gfx modify g_element "/" point glyph axes_xyz general size "50*50*50" select_on material blue selected_material default

# open Graphics Window
gfx cre win

# open Scene Editor
gfx edit scene
