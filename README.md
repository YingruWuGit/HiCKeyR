# HiCKeyR

R package of HiCKey. It has two functions ```segment(argv)``` and ```segHeatMap(argv, s, e)```.

# Install

First, user needs Rcpp and devtools:
```
install.packages("Rcpp")
install.packages("devtools")
library(devtools)
```
Then, install the HiCKeyR package:
```
install_github("https://github.com/YingruWuGit/HiCKeyR.git")
library(HiCKeyR)
```

# segment(argv)

This function is exactly the same as executing the C++ program. It does not return a value, but generates an output file containing boundaries (change-points), hierarchical orders and p-values. If the HiC data file has name "xxxx", then the output file is named as "xxxx_output" in the same directory.

- ```argv``` This parameter is the full path and name of "arguments_HiCKey.txt" containing the arguments setting. Please refer to the HiCKey repository (https://github.com/YingruWuGit/HiCKey) for the arguments setting in "arguments_HiCKey.txt".

To use it:

Download "BrownianP.txt" from 

```
segment("C:/Users/Andrew/Documents/GitHub/HiCKeyR/arguments_HiCKey.txt")
```

# segHeatMap(argv, s, e)

This function is almost the same as the function ```segment(argv)```. The only difference is that it returns a submatrix of the HiC data in which the upper triangular part is the original HiC data but the lower triangular part is constructed from the detected change-points and their hierarchical orders. User can draw a heatmap by the return matrix.

- ```argv``` Same as above.
- ```s``` Start index of the submatrix. Its default value is 0. Note that the index is 0 based.
- ```e``` End index of the submatrix. Its default value is -1, which means the end of the whole HiC matrix. Note that the last index of the submatrix is actually e-1 following C++ convention.
- ```return``` It returns a Rcpp::NumericMatrix from s to e. Its upper triangular part is the original HiC data but the lower triangular part is constructed from the detected change-points and their hierarchical orders. It also generates an output file containing change-point locations, orders and p-values.

For example:

If the HiC data is "samp_nested.txt" (https://github.com/YingruWuGit/HiCKey/tree/master/examples), then the "arguments_HiCKey.txt" can be
```
C:/Users/Andrew/Documents/GitHub/HiCKey/examples/samp_nested.txt
C:/Users/Andrew/Documents/GitHub/HiCKey/BrownianP.txt
m
3
0.05
0.00005
```
Input the following and leave s and e as default
```
X = segHeatMap("C:/Users/Andrew/Documents/GitHub/HiCKeyR/arguments_HiCKey.txt")
heatmap(X, scale = "none", Rowv = NA, Colv = NA, col = gray.colors(50, start = 1, end = 0, gamma = 0.15))
```
The heatmap is

![alt text](https://github.com/YingruWuGit/HiCKeyR/blob/main/sample_heatmap.png)

If the HiC data is "chr21_50kb.RAWobserved" (https://github.com/YingruWuGit/HiCKey/tree/master/examples) with resolution 50k, then input something like:
```
X = segHeatMap("C:/Users/Andrew/Documents/GitHub/HiCKey/arguments_HiCKey.txt", 0, 15050000)
heatmap(X, scale = "none", Rowv = NA, Colv = NA, col = gray.colors(50, start = 1, end = 0, gamma = 0.15))
```
