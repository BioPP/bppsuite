%define _basename bppsuite
%define _version 2.2.0
%define _release 1
%define _prefix /usr

URL: http://home.gna.org/bppsuite/

Name: %{_basename}
Version: %{_version}
Release: %{_release}
License: CECILL-2.0
Vendor: The Bio++ Project
Source: http://biopp.univ-montp2.fr/repos/sources/%{_basename}-%{_version}.tar.gz
Summary: The Bio++ Program Suite
Group: Productivity/Scientific/Other

Requires: libbpp-phyl9 = %{_version}
Requires: libbpp-seq9 = %{_version}
Requires: libbpp-core2 = %{_version}

BuildRoot: %{_builddir}/%{_basename}-root
BuildRequires: cmake >= 2.6.0
BuildRequires: gcc-c++ >= 4.0.0
BuildRequires: groff
BuildRequires: texinfo >= 4.0.0
BuildRequires: libbpp-core2 = %{_version}
BuildRequires: libbpp-core-devel = %{_version}
BuildRequires: libbpp-seq9 = %{_version}
BuildRequires: libbpp-seq-devel = %{_version}
BuildRequires: libbpp-phyl9 = %{_version}
BuildRequires: libbpp-phyl-devel = %{_version}


AutoReq: yes
AutoProv: yes
%if 0%{?mdkversion} >= 201100 || %{?distribution} == "Mageia"
BuildRequires: xz
%define zipext xz
%else
%if 0%{?mdkversion}
BuildRequires: lzma
%define zipext lzma
%else
BuildRequires: gzip
%define zipext gz
%endif
%endif

%description
Bio++ program suite includes programs:
 - BppML for maximum likelihood analysis,
 - BppSeqGen for sequences simulation,
 - BppAncestor for ancestral states reconstruction,
 - BppDist for distance methods,
 - BppPars for parsimony analysis,
 - BppSeqMan for file conversion and sequence manipulation,
 - BppConsense for building consensus tree and computing bootstrap values,
 - BppPhySamp for phylogenetic sampling,
 - BppReRoot for tree rerooting.
 - BppTreeDraw for tree drawing.
 - BppAlnScore for comparing alignments and computing alignment scores.
 - BppMixedLikelioods for computing the site per site likelihoods of submodels from a mixture model.
 
%prep
%setup -q

%build
CFLAGS="-I%{_prefix}/include $RPM_OPT_FLAGS"
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=%{_prefix}"
if [ %{_lib} == 'lib64' ] ; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DLIB_SUFFIX=64"
fi
if [ %{zipext} == 'lzma' ] ; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DDOC_COMPRESS=lzma -DDOC_COMPRESS_EXT=lzma"
fi
if [ %{zipext} == 'xz' ] ; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DDOC_COMPRESS=xz -DDOC_COMPRESS_EXT=xz"
fi

cmake $CMAKE_FLAGS .
make
make info

%install
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig

%files
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%{_prefix}/bin/bppml
%{_prefix}/bin/bppseqgen
%{_prefix}/bin/bppancestor
%{_prefix}/bin/bppdist
%{_prefix}/bin/bpppars
%{_prefix}/bin/bppseqman
%{_prefix}/bin/bppconsense
%{_prefix}/bin/bppphysamp
%{_prefix}/bin/bppreroot
%{_prefix}/bin/bpptreedraw
%{_prefix}/bin/bppalnscore
%{_prefix}/bin/bppmixedlikelihoods
%{_prefix}/share/info/bppsuite.info.%{zipext}
%{_prefix}/share/man/man1/bppml.1.%{zipext}
%{_prefix}/share/man/man1/bppseqgen.1.%{zipext}
%{_prefix}/share/man/man1/bppancestor.1.%{zipext}
%{_prefix}/share/man/man1/bpppars.1.%{zipext}
%{_prefix}/share/man/man1/bppdist.1.%{zipext}
%{_prefix}/share/man/man1/bppconsense.1.%{zipext}
%{_prefix}/share/man/man1/bppseqman.1.%{zipext}
%{_prefix}/share/man/man1/bppreroot.1.%{zipext}
%{_prefix}/share/man/man1/bppphysamp.1.%{zipext}
%{_prefix}/share/man/man1/bpptreedraw.1.%{zipext}
%{_prefix}/share/man/man1/bppalnscore.1.%{zipext}
%{_prefix}/share/man/man1/bppmixedlikelihoods.1.%{zipext}

%changelog
* Mon Sep 28 2014 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.2.0-1
- Compatibility update. Bio++ Program Suite version number is now indexed
  on Bio++'s version.
- Programs support the --seed argument for setting the random seed.
- bppSeqGen suport generic characters as input.
- bppPhySamp outputs sampled trees.
* Fri Mar 08 2013 Julien Dutheil <julien.dutheil@univ-montp2.fr> 0.8.0-1
- New models for proteins (COaLA)
- New program bppMixedLikelihoods
* Wed Feb 15 2012 Julien Dutheil <julien.dutheil@univ-montp2.fr> 0.7.0-1
- More models, sequence formats and bugs fixed. New bppAlnScore program.
* Thu Jun 09 2011 Julien Dutheil <julien.dutheil@univ-montp2.fr> 0.6.2-1
* Mon Feb 28 2011 Julien Dutheil <julien.dutheil@univ-montp2.fr> 0.6.1-1
* Mon Feb 07 2011 Julien Dutheil <julien.dutheil@univ-montp2.fr> 0.6.0-1
* Thu Mar 25 2010 Julien Dutheil <julien.dutheil@univ-montp2.fr> 0.5.0-1
* Wed Jun 10 2009 Julien Dutheil <jdutheil@birc.au.dk> 0.4.0-1
* Thu Dec 11 2008 Julien Dutheil <jdutheil@birc.au.dk> 0.3.1-1
* Thu Sep 23 2008 Julien Dutheil <jdutheil@birc.au.dk> 0.3.0-1
- Initial spec file.
