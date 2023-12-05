# BppSuite presentation

BppSuite is a suite of ready-to-use programs for phylogenetic and sequence analysis. 

## Installation 

### Standalone executables 

Standalone executables are available for [linux64](https://github.com/BioPP/bppsuite/releases/tag/v2.3.2)

[//]: [win32](http://biopp.univ-montp2.fr/repos/exe/win32/), [win64](http://biopp.univ-montp2.fr/repos/exe/win64/) and [Mac](http://biopp.univ-montp2.fr/repos/exe/mac/)

### From source files 

#### Get the sources 

This is done with <tt>git</tt>, for example in directory <tt>$bpp_dir</tt>:

<h5>
<pre>
cd $bpp_dir
git clone https://github.com/BioPP/bppsuite
</pre>
</h5>

#### Compiling 

Bio++ libraries need to be installed beforehand, for example in <tt>$bpp_dir</tt>. The needed libraries are [bpp-core](https://github.com/BioPP/bpp-core), [bpp-seq](https://github.com/BioPP/bpp-seq), [bpp-phyl](https://github.com/BioPP/bpp-phyl), [bpp-popgen](https://github.com/BioPP/bpp-popgen).

After, you proceed:

<h5>
<pre>
cd bppsuite
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir ./ # prepare compilation
make # compile
make install # move files to the installation directory (this will create a $bpp_dir/bin/ directory)
</pre>
</h5>

That's it ! The executables are now installed in <tt>$bpp_dir/bin</tt>. 

Without the option <tt>-DCMAKE_INSTALL_PREFIX=$bpp_dir</tt>, the standard <tt>/usr/local</tt> directory will be used, and the executables installed in  <tt>/usr/local/bin</tt>, a location which requires superuser access rights.

## Usage

Bppsuite executables should know where the dynamic libraries are.  A way to check it is the command:

<h5>
<pre>
ldd $bpp_dir$/bin/bppml
</pre>
</h5>

To configure this, set in the shell environment variable :

<h5>
<pre>
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$bpp_dir
</pre>
</h5>

(and source the configuration file or relog).

## Documentation

You can generate the pdf documentation by typing :

<h5>
<pre>
make pdf
</pre>
</h5>

and in html by typing:

<h5>
<pre>
make html
</pre>
</h5>

### Examples

Many examples are available in the subdirectory of <tt>Examples</tt>.

### Documentation 

Documentation can be found at: https://pbil.univ-lyon1.fr/bpp-doc/bppsuite/bppsuite.html

 
