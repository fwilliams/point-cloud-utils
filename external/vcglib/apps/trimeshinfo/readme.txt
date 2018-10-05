
   VCGLib  http://www.vcglib.net                                         o o     
   Visual and Computer Graphics Library                            o     o   
                                                                  _   O  _   
   Copyright(C) 2004-2006                                           \/)\/    
   Visual Computing Lab  http://vcg.isti.cnr.it                    /\/|      
   ISTI - Italian National Research Council                           |      
                                                                      \      
   TriMeshInfo 1.22 2006/10/30
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

* Number of Vertices 
* Number of Faces 
* Number of Edges 
* Number of Connected Components 
* Number of Boundaries 
* Number of Isolated Vertices (i.e. Unreferenced)
* Number of Duplicated vertices (duplicated vertices are referenced vertices
                                 which have the same positon in the space)
* Manifold 
* Genus (computed only for Manifold Datasets) 
* Self-Intersection (currently computed only for Datasets with less than 3M faces) 
* Orientability 
* Orientation
* Regularity (We consider REGULAR those meshes that have 6 incident edges
  for each internal vertex, and 4 incident edges for each vertex on the 
  boundary. In all other cases we consider the mesh irregular.)
* Volume
* Attribute present on mesh (colors, texture coords etc);

The application has no graphical interface but works as the "Metro" tool on command line. 

TriMeshInfo is written in C++ and makes use of the VCG library. 
The tool supports the following file formats:

  - PLY 
  - OFF (http://www.geomview.org/docs/html/geomview_41.html#SEC44)
  - STL (http://astronomy.swin.edu.au/~pbourke/dataformats/stl/) 
  - OBJ (http://www.eg-models.de/formats/Format_Obj.html)


--- Command-line Reference ---

Usage:

    TriMeshInfo <mesh> [options]

Valid options are the following:

    -q Quiet (disable verbose mode that is enabled by default)
    -x Enable XML output
    -h Enable HTML output
    -s <filename> Save the clean mesh

The HTML output creates in the directory where TriMeshInfo is launched a file
called "result.html". This file contains an hmtl table with the retrieved 
mesh information. 
If this file is just present in the working directory the output of the TriMeshInfo 
is added to the existing table. In this way it is possible to summarize the results
obtained from several meshes.

If you choose to save the "clean" mesh, the mesh without its unreferenced vertices
and with the duplicated vertices merged is saved.

