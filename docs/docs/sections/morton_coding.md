# Morton Encoding/Decoding a Point Cloud
Point Cloud Utils has a number of tools for operating on Morton Encoded point clouds. Morton encoding works by projecting points onto a space filling [Z-Curve](https://en.wikipedia.org/wiki/Z-order_curve) (illustrated below), and recording the distance along this curve from the origin. 

<p align="center">
  <div class="row" style='content: "";clear: both; display: table;'>
    <div class="column" style="float: left; width: 33.33%; padding: 5px;">
      <img src="../../imgs/z_curve_2.png" alt="Snow" style="width:100%">
    </div>
    <div class="column" style="float: left; width: 33.33%; padding: 5px;">
      <img src="../../imgs/z_curve_4.png" alt="Forest" style="width:100%">
    </div>
    <div class="column" style="float: left; width: 33.33%; padding: 5px;">
      <img src="../../imgs/z_curve_8.png" alt="Mountains" style="width:100%">
    </div>
  </div> 
  <figcaption style="text-align: center; font-style: italic;">Z-order curves fill space by repeated subdivision. This image shows one, two and three iterations of z-order subdivision.</figcaption> 
</p>

Morton Encoding point clouds has a number of useful application, such as approximate-k-nearest neighbor search, point hashing, and point sorting, to name a few.

## Morton Encoding and Decoding Points

## Approximate K-Nearest-Neighbor Search with Morton Coding

## Sorting Points Along a Morton Curve

## Translating Points in Morton Space

