
   VCGLib  http://www.vcglib.net                                         o o     
   Visual and Computer Graphics Library                            o     o   
                                                                  _   O  _   
   Copyright(C) 2005-2006                                           \/)\/    
   Visual Computing Lab  http://vcg.isti.cnr.it                    /\/|      
   ISTI - Italian National Research Council                           |      
                                                                      \      
   Metro 4.07 2007/05/11 
   All rights reserved.                                                      
   
                                                                       
This program is free software; you can redistribute it and/or modify      
it under the terms of the GNU General Public License as published by      
the Free Software Foundation; either version 2 of the License, or         
(at your option) any later version.                                       
                                                                          
This program is distributed in the hope that it will be useful,           
but WITHOUT ANY WARRANTY; without even the implied warranty of            
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             
GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          
for more details.                                                 

--- Synopsis ---

Metro is a tool designed to evaluate the difference between two triangular meshes. 
Metro adopts an approximated approach based on surface sampling and point-to-surface distance computation. 
Please, when using this tool, cite the following reference:


P. Cignoni, C. Rocchini and R. Scopigno
"Metro: measuring error on simplified surfaces"
Computer Graphics Forum, Blackwell Publishers, vol. 17(2), June 1998, pp 167-174


For any question about this software please contact:
Paolo Cignoni ( paolo.cignoni@isti.cnr.it )

--- General Info ---

Metro is a tool designed to evaluate the difference between two triangular meshes. 
Metro adopts an approximated approach based on surface sampling and point-to-surface distance computation. 
Three different surface sampling methods are implemented:

    *   Montecarlo sampling (pick k random samples in the interior of each face)
    *   Subdivision sampling (recursively subdivide each face along the longest edge and choose the sample in the center of each cell)
    *   Similar Triangles sampling (subdivide each face F in k polygons similar to F and sample the face in correspondence with the vertices of these polygons, internal to F)

Note that the three methods described above are used to sample only the interior of each face. 
A different scheme is used to sample vertices and edges: vertices are sampled in the straightforward manner, 
while edges are sampled by uniformly interleaving samples along each edge.

Three different Spatial indexing structures can be used to find the closest point to a sample, a Statically Allocated Uniform Grid, a Hashed Uniform Grid and a Hierarchy of axis aligned bounding boxes.

--- Basic usage ---

Metro is a command-line tool which allows the user to select among different sampling schemes. 
A list of the command-line parameters accepted by the tool is shown in the following.

Usage: Metro file1 file2 [opts]

where "file1" and "file2" are the input meshes in PLY, OFF or STL format, and opts can be:

  -v         disable vertex sampling
  -e         disable edge sampling
  -f         disable face sampling
  -u         ignore unreferred vertices
  -sx        set the face sampling mode
             where x can be:
              -s0  montecarlo sampling
              -s1  subdivision sampling
              -s2  similar triangles sampling (Default)
  -n#        set the required number of samples (overrides -A)
  -a#        set the required number of samples per area unit (overrides -N)
  -c         save a mesh with error as per-vertex colour and quality
  -C # #     Set the min/max values used for color mapping
  -L         Remove duplicated and unreferenced vertices before processing
  -h         write files with histograms of error distribution
  -G         Use a static Uniform Grid as Search Structure (default)
  -A         Use an Axis Aligned Bounding Box Tree as Search Structure
  -H         Use an Hashed Uniform Grid as Search Structure
  -O         Use an Octree as Search Structure
  
  
The -C option is useful in combination with -c option for creating a set of 
meshes with a coherent coloring scheme. 
It sets how the errors are mapped into color according to the following formula, 
let e be the error and ColorRamp be a R->RGB function mapping 0..1 values 
into a smooth RedYellowGreenCyanBlue ramp:

			e=Clamp(e,min,max);
			VertexColor = ColorRamp( (e-min)/(max-min) );

The Histogram files saved by the -h option contains two column of numbers 
e_i and p_i; p_i denotes the fraction of the surface having an error 
between e_i and e_{i+1}. The sum of the second column values should give 1.