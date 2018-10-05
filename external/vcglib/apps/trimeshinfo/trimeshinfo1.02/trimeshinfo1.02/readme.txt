
   VCGLib  http://www.vcglib.net                                         o o     
   Visual and Computer Graphics Library                            o     o   
                                                                  _   O  _   
   Copyright(C) 2004-2005                                           \/)\/    
   Visual Computing Lab  http://vcg.isti.cnr.it                    /\/|      
   ISTI - Italian National Research Council                           |      
                                                                      \      
   TriMeshInfo 1.0 21/09/2004
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

TriMeshInfo is a tool designed to inspect 3D models and retrieve all the 
topological related information. It can be used to automate the process 
of decoding 3D mesh inherent properties and ease data classification 
and retrieval. 


For each analyzed dataset the following information are extracted: 

* Number of Vertices (Unreferenced vertices are listed separately) 
* Number of Faces 
* Number of Edges 
* Number of Connected Components 
* Number of Boundaries 
* Number of Isolated Vertices (i.e. Unreferenced) 
* Manifold 
* Genus (Computed only for Manifold Datasets) 
* Self-Intersection (Currently computed only for Datasets with less than 3M faces) 
* Orientability 
* Orientation 
* Regularity (We consider as regular those meshes generated through 
  regular subdivision. Each non boundary vertex of a regular mesh has 
  6 incident edges, if there are only 5 incident edges the mesh is said to be 
  semi-regular, irregular in all other cases) 

The application has no graphical interface but works as the "Metro" tool on command line. 

TriMeshInfo is written in C++ and makes use of the VCL library. 
The tool supports two file formats ply and off (as described in http://www.geomview.org/docs/html/geomview_41.html#SEC44) . 
