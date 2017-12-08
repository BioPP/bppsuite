[[Category:BppSuite]]

BppSuite is a suite of ready-to-use programs for phylogenetic and sequence analysis. 

[http://biopp.univ-montp2.fr/forge/bppsuite Here] is the forge where the programs are described.

== Installation ==

=== Standalone executables ===

Standalone executables are available for [http://biopp.univ-montp2.fr/repos/exe/lin32/ linux32], [http://biopp.univ-montp2.fr/repos/exe/lin64/ linux64], [http://biopp.univ-montp2.fr/repos/exe/win32/ win32], [http://biopp.univ-montp2.fr/repos/exe/win64/ win64] and [http://biopp.univ-montp2.fr/repos/exe/mac/ Mac].

=== From source files ===

===== Getting the sources =====

This is done with <tt>git</tt>, for example in directory <tt>$bpp_dir</tt>:

<pre>
cd $bpp_dir
git clone https://github.com/BioPP/bppsuite
</pre>

===== Compiling =====

Bio++ libraries shoud be installed beforehand, for example in <tt>$bpp_dir</tt>.

<pre>
cd bppsuite
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir ./ # prepare compilation
make # compile
make install # move files to the installation directory (this will create a $bpp_dir/bin/ directory)
</pre>

That's it ! The executables are now installed in <tt>$bpp_dir/bin</tt>. 
For more information on how to compile and run Bio++ dependent programs, see the [[:Category:Usage|Usage]] pages.

* Without the option <tt>-DCMAKE_INSTALL_PREFIX=$bpp_dir</tt>, the standard <tt>/usr/local</tt> directory will be used, and the executables installed in  <tt>/usr/local/bin</tt>, a location which requires superuser access rights.

===== Usage =====

Bppsuite executables should know where the dynamic libraries are.  A way to check it is the command:

<pre>
ldd $bpp_dir$/bin/bppml
</pre>

To configure this, set in the .bashrc the environment variable :

<pre>
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$bpp_dir
</pre>

(and source this  file or relog).

===== Documentation =====

You can also generate the pdf documentation by typing :

<pre>
make pdf
</pre>

== Examples ==

Many examples are available in the subdirectory of <tt>Examples</tt>.

== Documentation ==

Documentation can be found at [http://biopp.univ-montp2.fr/manual/html/bppsuite/ http://biopp.univ-montp2.fr/manual/html/bppsuite/].
 