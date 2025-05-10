# Asymptotic and Regime Regression Program (ARRP)
**ARRP** implements both regime regression and robust regression corresponding to linear and constant models. 
* **Robust regression** is a systematic method of removing outliers from linear regression data.
* **Regime regression** is a systematic method for removing points from regression data when they lie outside the contiguous regime in which a linear model is valid.
  * **Asymptotic regression** is a special type of regime regression, used for estimating parameters in an asymptotic regime.  
  * Generally, 'regime regression' is a better name than 'asymptotic regression' for the purpose of the ARRP.

---
### Files and Installation

The download ZIP file `ARRP_1.1.zip` contains a single directory `ARRP_1.1` containing:
1. `README.md`: contains this README.
1. `arrp_cpp_files/`: directory with C++ source files. `sls_arrp.t.cpp` contains main().
1. `arrp_1.1_LINUX.zip`: contains a single LINUX executable `arrp.exe` with the regression methods.
1. `arrp_1.1_WINDOWS.zip`: contains a single Windows executable `arrp.exe` with the regression methods.
1. `points.dat`: an example of an input file.

* No special installation is required.
* The executable files can be downloaded, unzipped and run with the appropriate command line.
* Alternatively, the source C++ files can be downloaded and complied in a suitable C++ environment.

---
### Usage

**`arrp.exe <OPTIONS>`**

* `arrp.exe` is the name of the binary file.
* `<OPTIONS>` is the list of options.
  * Options requiring a parameter value have the form  
`—<option> <parameter value>`.
  * If the parameter value has a default, the option may be omitted, as indicated by square brackets below.   
`[ —<option> <parameter value> ]`.
  * Options not requiring a parameter value always have a default, so they have the form  
`[ —<option> ]`.
---
### Options

**Examples**

* `arrp.exe -in points.dat -out regress.out`  
Linear regime regression (c = 2.0), allowing removal of points on both left and right
* `arrp.exe -in points.dat -out regress.out -flat -include left`  
Constant regime regression (c = 2.0), onto a horizontal line, allowing removal of points on only the right
* `arrp.exe -in points.dat -out regress.out -c 1.0 -include right`  
Linear regime regression (c = 1.0), allowing removal of points on only the left
* `arrp.exe -in points.dat -out regress.out -regression weighted`  
Weighted linear regression
* `arrp.exe -in points.dat -out regress.out -regression robust`  
Linear robust regression with the Andrews' function (a = 1.339) [2]
* `arrp.exe -in points.dat -out regress.out -flat -regression robust -a 3.0`  
Constant robust regression with the Andrews' function (a = 3.0) [2]
---
### Options and Parameter Values

* `-in <Input file name>`  
The file contains the data points for the regression.
The file format is described below.
* `-out <Output file name>`  
The output file format is described below.
* `[ -append ]`  
Default: If the output file already exists, overwrite it with the new output.
  * `-append`: If the output file already exists, append output to it.   
The option is useful for collecting output from multiple related runs.
* `[ -flat ]`, an option for switching to constant ("flat") regression on a line of 0 slope.  
  * Default: If the option is omitted, ARRP defaults to linear regression y = β<sub>0</sub>+β<sub>1</sub>x.
* `-flat`: ARRP performs constant regression y = β<sub>0</sub>.
* `[ -regression <Type> ]`, an option for switching to robust or weighted least-square regression  
Default `-regression asymptotic`: If the option is omitted, ARRP defaults to asymptotic (regime) regression. 
ARRP then uses the following additional parameters for its regime regression:  
  * `[ -c <Parameter 'c'> ]`, the parameter 'c' of regime regression  
    * Default `'-c 2.0'`, as recommended in [1]  
    * The parameter 'c' must be strictly positive.
  * `[ -include <parameter> ]`  
    * Default `-include neither`: If the option is omitted, the regime regression permits removal of data points on both the left and right sides, not forcing the inclusion of points on either boundary.
    * `-include left`: the regime regression must include points on the left starting from the leftmost point, but permits removal of data points from the right side.
    * `-include right`: the regime regression includes points on the right starting from the rightmost point, but permits removal of data points from the left side.
* `-regression robust: ARRP performs robust regression.  
ARRP then uses the following additional parameters for its robust regression.
  * `[ -a <Parameter 'a'> ]`, the parameter 'a' of robust regression
  * Default `-a 1.339`, as recommended in [2]
  * The parameter `'a'` must be strictly positive.
  * Implementation notes:  
ARRP applies robust regression twice during its run:
      * to provide an initial approximation with Ramsay’s function (with `‘a’` = 0.3),
      * then, to do a robust regression with Andrews' function (with the actual input parameter `‘a’`).
* `-regression weighted`: the regime regression includes all points, making the regime regression equivalent to a weighted least-squares regression.
---
### The Format of an Input File
<pre>
x<sub>1</sub>    y<sub>1</sub>    error<sub>1</sub> 
x<sub>2</sub>    y<sub>2</sub>    error<sub>2</sub> 
................. 
x<sub>n</sub>    y<sub>n</sub>    error<sub>n</sub> 
</pre>

* <code>(x<sub>k</sub>, y<sub>k</sub>, error<sub>k</sub>)</code> (k = 1,2,...,n) are the data points (with errors) for regression.
* <code>error<sub>k</sub></code> is the error in <code>y<sub>k</sub></code>, so error<sub>k</sub> must be strictly positive.
* The points must be sorted, so <code>x<sub>1</sub> < x<sub>2</sub> <...< x<sub>n</sub></code>.
* The output file contains all ARRP's input parameters, including the input data filename. The input data can be annotated elsewhere if necessary, so the ARRP input data format is strict.
  * Empty lines and lines containing only whitespace are ignored.
  * Other lines contain exactly 1 point, i.e., 3 numbers, separated by 2 strings of whitespace.
  * All deviations from the input format are reported.
 
**The examples below use the following input file, "points.dat".**
<pre>
0       0.33726     0.0165   
0.1     0.38115     0.0165   
0.2     0.10626     0.0165   
0.3     1.056       0.0165   
0.4     1.03026     0.0165   
0.5     1.02696     0.0165   
0.6     1.0593      0.0165   
0.7     1.0527      0.0165   
0.8     1.03323     0.0165   
0.9     1.05336     0.0165 
</pre>

### The Format of an Output File

In the example output files below, all whitespace is either tab or newline.  
In them, all lines (including the last) end in a newline.

**An Example of the Output for Regime Regression**

`arrp.exe -in points.dat -out regress.out -flat -include right`  
The command does constant regime regression (c = 2.0), onto a horizontal line, allowing removal of points on the left only.

**"regress.out"**
<pre>
Input      points.dat 
Flat       yes        Result     f(y)~beta0
Regression asymptotic c          2          Include    right   
Range      0.3        0.9
beta0	   1.04454    error      0.00623641   
	
X          Y          Error      Weight
0          0.33726    0.0165     0
0.1        0.38115    0.0165     0
0.2        0.10626    0.0165     0
0.3        1.056      0.0165     1
0.4        1.03026    0.0165     1
0.5        1.02696    0.0165     1
0.6        1.0593     0.0165     1
0.7        1.0527     0.0165     1
0.8        1.03323    0.0165     1
0.9        1.05336    0.0165     1
</pre>

* The weight "1" indicates that the point was used in the regime regression; "0", not.
* See the next example for output without the `-flat` option.
* 
**An Example of the Output for Robust Regression**

`arrp.exe -in points.dat -out regress.out -regression robust`  
The command does robust linear regression (a = 1.339).

**"regress.out"**
<pre>
Input      points.dat 
Flat       no         Result     f(y)~beta0+beta1*y
Regression robust     a          1.339                     
Range      0          0.9
beta0	   0.320971   error      0.0123775  
beta1	   0.831244   error      0.0211017   
	
X          Y          Error      Weight
0          0.33726    0.0165     3349.24
0.1        0.38115    0.0165     3047.52
0.2        0.10626    0.0165     0
0.3        1.056      0.0165     0
0.4        1.03026    0.0165     0
0.5        1.02696    0.0165     0
0.6        1.0593     0.0165     0
0.7        1.0527     0.0165     0
0.8        1.03323    0.0165     1446.95
0.9        1.05336    0.0165     3370.52
</pre>

* A data point with 0 weight is considered an outlier and does not influence the regression result.  
  
* The -append option prefixes a delimiter to the output above.
* To demarcate appended output, the -append option prefixes a line of hyphens.

---
### References

1. S.L. Sheetlin,Y. Park",J.L. Spouge (2011)  
An objective method for estimating asymptotic parameters, with application to biological sequence alignment  
Physics Review E 84 : 031914 https://link.aps.org/doi/10.1103/PhysRevE.84.031914
1. D.C. Montgomery,E.A. Peck, G.G. Vining (2001)  
Introduction to linear regression analysis
3rd edition, Wiley-Interscience

  
