We have the original geometry file named "Doua.geo_unrolled". The dimensions of this area are approximately 170m x 200m x 19m.

From this original geometry, we have created four physical surfaces:
- "Sol", representing the ground,
- "Ext", representing the external boundary,
- "Batiment", representing the buildings,
- "Gamma", regrouping "Sol", "Ext" and "Batiment".
Additionally, we have defined a physical volume named "Omega", which encompasses the entire area of interest.

These modifications have been applied in the file "Doua_phy_names.geo_unrolled".

Building upon the "Doua_phy_names.geo_unrolled" file, we increased the height between the ceiling and the ground. To achieve this, we modified the z-coordinates of four specific points from 20 to 50. These points are identified by their tags in the list (35, 36, 122, 123).

The alterations are recorded in the file "Doua_phy_names_hplus.geo_unrolled".