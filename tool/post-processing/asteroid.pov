// Asteroid
#include "colors.inc"
#include "shapes.inc"
#include "textures.inc"
#include "stones.inc"
#include "woods.inc"
#include "glass.inc"
#include "metals.inc"  
#include "functions.inc"
#version 3.6;
// ==============================================================
#declare BASE_SHAPE=
function
{
  sqrt(x*x+y*y+z*z) - 1
}

#declare CRATER_SHAPE_TEMPLT=
function
{
  pigment
  {
    crackle form <1.5,0,0>
    color_map
    {
      [0 rgb <1.0,1.0,1.0>]
      [0.75 rgb <0.0,0.0,0.0>]
      [1 rgb <0.2,0.2,0.2>]
    }
    cubic_wave
  }
}

#declare CRATER_SHAPE=
function(x,y,z,S)
{
  CRATER_SHAPE_TEMPLT(x/S,y/S,z/S).red
}

#declare asteroid_regolith=
  texture   { pigment {  bozo    color_map { [0 rgb <0.3,0.3,0.3>] [1 rgb <1.0,1.0,1.0>] }  scale 0.2 
                         warp { turbulence .5  octaves 3  omega 1.0 lambda .7 } scale 0.5 }
                 finish {ambient 0.0 diffuse 1.0 brilliance 1.0 specular 0.1 roughness 0.08 } 
                 normal { agate 0.13 scale 0.08 }   }
                 
/*
isosurface{
  function{
      BASE_SHAPE(x, y, z)    
      + f_noise3d(x,y,z)
      + .04 * CRATER_SHAPE(x,y,z,.35) 
      + .015 * CRATER_SHAPE(x+10,y+10,z+10,.15)
      + .015 * CRATER_SHAPE(x,y,z,.1)
      + .005 * CRATER_SHAPE(z+1,x+3,y+2,.05)

  }
  contained_by{box{-1.2 1.2}}
  threshold 0
  texture{asteroid_regolith}
}
*/

sphere {<   0.0,  0.0,   0.0>,    0.5 
    texture {asteroid_regolith}
}   
        
camera {location <  0.0000e+00,  -5.0000e+00,  0.0000e+0> 
        sky   <  0.0000e+00,  0.0000e+00,  1.0000e+00> 
        look_at <  0.0000e+0,  0.0000e+00,  0.0000e+0> 
        up <  0.0000e+00,  9.0000e-01,  0.0000e+00> 
        right <  1.6000e+00,  0.0000e+00,  0.0000e+00> 
        }   
light_source {<  1.5000e+00,  0.0000e+00,  1.5000e+01> color 1.5}
light_source {<  0.0000e+00,  -5.0000e+00,  0.0000e+0> color 0.5}  
