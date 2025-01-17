# CS 73/273 Final Project: Gradient Domain Copy-Paste
## Mark Gitau, Sihau Huang, Ryan Tucker, and Irene Zuo

This project uses base code from Dartmouth's Fall 2020 offering of CS 73/273 - "Computational Aspects of Digital Photography".
 
See course [website](https://canvas.dartmouth.edu/courses/43075) for details.

Final presentation [website](https://docs.google.com/presentation/d/1mRlbEP9puaESHQSG6Zz9kIsjdmanST-bD5DV9bQyqQM/edit?usp=sharing)

## Credits
The basecode is written by Wojciech Jarosz, but it is heavily derived from (with permission):

:information_source: MIT's 6.815/6.865 basecode, written and designed by:
* Frédo Durand
* Katherine L. Bouman
* Gaurav Chaurasia
* Adrian Vasile Dalca
* Neal Wadhwa

### References
* Tutorial by Eric Arnebäck https://erkaman.github.io/posts/poisson_blending.html
* Poisson Image Editing, Patrick Pe ́rez, Michel Gangnet and Andrew Blake, Microsoft Research UK

## Content
- All functions demonstrating Gradient Domain copy-paste are cointained within the file `a6-main.cpp`, namely `testBlend()`, `poisson_blending_gradient_descent()`, `testMixedBlend()`, `testFlatten()`, and `testColor()`.
1. `testBlend()`: a function which demonstrates basic gradient domain copy-pasting using Poisson cloning.
2. `poisson_blending_gradient_descent()`: a function which demonstrates gradient domain copy-pasting using gradient descent, operated in the log domain.
3. `testMixedBlend()`: a function which demonstrates mixed blending, where the gradient for each pixel is calculated using both source and target image gradients.
4. `testFlatten()`: a function which demonstrates flattening, where the gradient for each pixel is only calculated across edges, giving the output a 'flat' appearance.
5. `testColor()`: a function which demonstrates local color changing, where some part of the input image's color is altered, and the part is blended back into the image.
6. `testLumi()`: a function which demonstrates local luminance changing.
- The code makes use of the Eigen library for vectors, matrices, and solvers.
