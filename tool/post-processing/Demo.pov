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
camera {location <  1.2500e+00,  5.0000e+00,  5.0000e+0> 
        sky   <  0.0000e+00,  0.0000e+00,  1.0000e+00> 
        look_at <  1.2500e+0,  0.0000e+00,  0.0000e+0> 
        up <  0.0000e+00,  9.0000e-01,  0.0000e+00> 
        right <  1.6000e+00,  0.0000e+00,  0.0000e+00> 
        }   
light_source {<  10.5000e+00,  20.0000e+00,  5.0000e+01> color 1.5}
light_source {<  0.0000e+00,  20.0000e+00,  0.0000e+0> color 0.2}  
// ==============================================================
// =======================Asteroid===============================
// ==============================================================
#declare BASE_SHAPE=
function
{
  sqrt(x*x+y*y+z*z) - 0.5
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

#declare asteroid=
  texture   { pigment {  bozo    color_map { [0 rgb <0.3,0.3,0.3>] [1 rgb <1.0,1.0,1.0>] }  scale 0.2 
                         warp { turbulence .5  octaves 3  omega 1.0 lambda .7 } scale 0.5 }
                 finish {ambient 0.0 diffuse 1.0 brilliance 1.0 specular 0.1 roughness 0.08 } 
                 normal { agate 0.13 scale 0.08 }   }
                 
isosurface{
  function{
      BASE_SHAPE(x, y, z)    
      + .04 * CRATER_SHAPE(x,y,z,.35) 
      + .015 * CRATER_SHAPE(x+10,y+10,z+10,.15)
      + .015 * CRATER_SHAPE(x,y,z,.1)
      + .005 * CRATER_SHAPE(z+1,x+3,y+2,.05)

  }
  contained_by{box{-1.2 1.2}}
  threshold 0
  texture{asteroid}
}
// ==============================================================
// =======================Asteroid Regolith======================
// ==============================================================
#declare asteroid_regolith_1= 
  texture   { pigment {  crackle form <1.0,0,0>  color_map { [0 rgb <0.15,0.15,0.15>] [1 rgb <0.3,0.3,0.3>] } cubic_wave  scale 0.1 
                         warp { turbulence .5  octaves 3  omega 1.0 lambda .7 } scale 0.5 } 
              finish {ambient 0.0 diffuse 1.0 brilliance 2.0 specular 0.05 roughness 0.1 } 
              normal { agate 0.13 scale 0.08 }   }
sphere {<   1.25,  0.0,   0.0>,    0.5 
    texture {asteroid_regolith_1}
}   
// ==============================================================
// =======================Asteroid Regolith======================
// ==============================================================
      
#declare White_Marble_Map = 
color_map { 
    [0.0 rgb <0.5, 0.5, 0.5>] 
    [0.8 rgb <0.3, 0.3, 0.3>] 
    [1.0 rgb <0.1, 0.1, 0.1>] 
} 
#declare White_Marble = 
pigment { 
    marble    turbulence 1.0    color_map { White_Marble_Map } 
} 
#declare asteroid_regolith_2= 
  texture {pigment { White_Marble} finish {ambient 0.0 diffuse 0.7 brilliance 2.0 specular 0.05 roughness 0.1} normal { agate 0.13 scale 0.08 }} 
#declare Rubber = 
  texture{ pigment{ aoi color_map{ 
	[0.00 rgb <.0075, .0175, .0025>] 
	[0.55 rgb <.020, .022, .024>] 
	[0.65 rgb <.004, .004, .004>] 
	[0.85 rgb <.006, .002, .001>] 
	[1.00 rgb <.007, .004, .001>]} 
    poly_wave 1.25  scale 0.1  } 
    normal {bozo 0.2 scale 0.25} 
    finish{ specular .015  roughness .075  brilliance 0.275 
    } 
}   
sphere {<   2.5,  0.0,   0.0>,    0.5 
    texture {asteroid_regolith_2 scale 0.2}}
// ==============================================================
// =======================Robot==================================
// ==============================================================
#declare PaintColor = color SkyBlue;
 
#declare PaintBright =
 
    pigment{
         PaintColor
        }
 
#declare PaintDark =
 
    pigment{
         PaintColor / 10
        }
 
#declare CarPaint = 
  texture{
    pigment{
         aoi
         pigment_map{
           //[0.0 PaintDark]
           [0.5 PaintDark]
           [1.0 PaintBright]
         }
        }
    normal {bozo 0.05 scale 0.1}
    finish{
         specular 0.4
         roughness 0.05 
         metallic
         diffuse 0.25
         brilliance 1
         reflection{ 0.05   
         }
    }
  }  
sphere {<   3.75,  0.0,   0.0>,    0.5 
    texture {CarPaint}}
    
// ==============================================================
// =======================Robot==================================
// ==============================================================
#declare PaintColor = color White;
 
#declare PaintBright =
 
    pigment{
         PaintColor
        }
 
#declare PaintDark =
 
    pigment{
         PaintColor / 2
        }
 
#declare CarPaint = 
  texture{
    pigment{
         aoi
         pigment_map{
           //[0.0 PaintDark]
           [0.5 PaintDark]
           [1.0 PaintBright]
         }
        }
    normal {bozo 0.05 scale 0.1}
    finish{
         diffuse 0.4
         brilliance 0.6
         reflection{
                rgb <.05, .05, .05>, rgb<.2,.2,.2>
                fresnel on       
         }
    }
  }
plane { <0, 0, 1>, -0.5
    texture {CarPaint}
}
// ==============================================================
// =======================Robot==================================
// ==============================================================
#declare Leather =
 
  texture{
    pigment{
      crackle
      metric 3
      turbulence .25
      color_map{
	[0.00 rgb <.015, .008, .004>]
	[0.25 rgb <.017, .008, .004>]
	[0.50 rgb <.010, .002, .004>]
	[0.75 rgb <.008, .002, .004>]
	[1.00 rgb <.006, .002, .003>]
      }
    scale 0.1
    }
    normal{
      crackle
      metric 3
      turbulence .25
      //granite .16
      scale 0.1
    }
    finish{
      specular .04
      roughness .1
    }
  }
sphere {<   5.0,  0.0,   0.0>,    0.5 
    texture {Leather}}
// ==============================================================
// =======================Sand===================================
// ==============================================================    
#declare Sand=
  texture{  pigment{ color rgb <.518, .339, .138> }
    normal{ bumps 5 scale 0.05 }
    finish{ specular .3 roughness .8 }
  }
 
  texture{ pigment{ wrinkles scale 0.05 color_map{
	[0.0 color rgbt <1, .847, .644, 0>]
	[0.2 color rgbt <.658, .456, .270, 1>]
	[0.4 color rgbt <.270, .191, .067, .25>]
	[0.6 color rgbt <.947, .723, .468, 0>]
	[0.8 color rgbt <.356, .250, .047, 1>]
	[1.0 color rgbt <.171, .136, .1, 1>]
      }
    }
  }
sphere {<   -1.25,  0.0,   0.0>,    0.5 
    texture {Sand scale 0.1}}
// ==============================================================
// =======================Snow===================================
// ==============================================================     
#declare Snow =
  texture{ 
    pigment{
      color rgb <0.9, 0.95, 1>
    }
    normal{
      gradient y .5
      slope_map {
	[0 <.1, .1>]
	[0.25 <.25, 0>] 
	[0.5 <.1, -.1>] 
	[0.75 <.25, 0>] 
	[1 <.1, .1>] 
      }
      turbulence 0.5
      scale 3
    }
    finish{
      brilliance 0.75
      phong 0.1
      phong_size 5
      subsurface {translucency <0.1, 0.31, 0.48>}
      //emission .2
      //use with  radiosity instead
    }
  }
 
  texture{
    pigment{
      color rgbt <1, 1, 1, .9>
    }
    normal{
      bumps 5
      scale .05  
    }
    finish{
      specular 1
      roughness .01
    }
  }
 
   texture{
    pigment{
      color rgbt <1, 1, 1, .9>
    }
    normal{
      bumps 3
      scale .1
    }
   }
sphere {<   -2.5,  0.0,   0.0>,    0.5 
    texture {Snow scale 0.2}}         
// ==============================================================
// =======================Stripe===================================
// ==============================================================     
#declare Stripe= 
  texture {pigment { 
               gradient <1,0,0> color_map {[0.0 color rgb<0,0,0>]
                                           [0.15 color rgb<0,0,0>]
                                           [0.15 color rgb<1,1,1>]
                                           [0.85 color rgb<1,1,1>]
                                           [0.85 color rgb<0,0,0>]
                                           [1.0 color rgb<0,0,0>]}
              scale 0.5}
finish {ambient 0.0 diffuse 0.7 brilliance 2.0 specular 0.05 roughness 0.1} normal { agate 0.13 scale 0.08 }}
sphere {<   1.25,  -1.5,   0.0>,    0.5 
    texture {Stripe scale 1.0}}