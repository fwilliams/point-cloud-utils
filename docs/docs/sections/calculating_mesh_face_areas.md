# Calculating Mesh Face Areas
You can easily calculate mesh face areas in Point-Cloud-Utils
```
import point_cloud_utils as pcu

v, f = pcu.load_mesh_vf("bunny.ply")

areas = pcu.mesh_face_areas(v, f)
```

<p align="center">
    <img src="../../imgs/face_areas.png" alt="Mesh colored by face areas" style="width:70%">
    <figcaption style="text-align: center; font-style: italic;">Mesh colored by face areas.</figcaption>
</p>