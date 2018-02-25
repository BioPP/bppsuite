%define _basename bppsuite
%define _version 2.4.0
%define _release 1
%define _prefix /usr

URL: https://github.com/BioPP

Name: %{_basename}
Version: %{_version}
Release: %{_release}
License: CECILL-2.0
Vendor: The Bio++ Project
Source: %{_basename}-%{_version}.tar.gz
Summary: The Bio++ Program Suite
Group: Productivity/Scientific/Other

Requires: libbpp-phyl11 = %{_version}
Requires: libbpp-seq11 = %{_version}
Requires: libbpp-core3 = %{_version}

BuildRoot: %{_builddir}/%{_basename}-root
BuildRequires: cmake >= 2.8.11
BuildRequires: gcc-c++ >= 4.7.0
BuildRequires: groff
BuildRequires: texinfo >= 4.0.0
BuildRequires: libbpp-core4 = %{_version}
BuildRequires: libbpp-core-devel = %{_version}
BuildRequires: libbpp-seq12 = %{_version}
BuildRequires: libbpp-seq-devel = %{_version}
BuildRequires: libbpp-phyl12 = %{_version}
BuildRequires: libbpp-phyl-devel = %{_version}
BuildRequires: libbpp-popgen7 = %{_version}
BuildRequires: libbpp-popgen-devel = %{_version}


AutoReq: yes
AutoProv: yes
%if 0%{?mdkversion}
%if 0%{?mdkversion} >= 201100
BuildRequires: xz
%define compress_program xz
%else
BuildRequires: lzma
%define compress_program lzma
%endif
%else
BuildRequires: gzip
%define compress_program gzip
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
 - BppReRoot for tree rerooting.
 - BppTreeDraw for tree drawing.
 - BppAlnScore for comparing alignments and computing alignment scores.
 - BppPopStats for population genetics.
 - BppMixedLikelioods for computing the site per site likelihoods of submodels from a mixture model.
 
%prep
%setup -q

%build
CFLAGS="$RPM_OPT_FLAGS"
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=%{_prefix} -DCOMPRESS_PROGRAM=%{compress_program}"
cmake $CMAKE_FLAGS .
make

%install
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig

%files
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%{_prefix}/bin/*
%{_prefix}/share/info/*.info*
%{_prefix}/share/man/man1/*.1*

%changelog
* Thu Feb 22 2018 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.4.0-1
* Tue Jun 06 2017 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.3.1-1
* Wed May 10 2017 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.3.0-1
- New BppPopStats program
- BppPhySamp is now distributed separately
- Several bugs fixed and improvements
* Mon Sep 29 2014 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.2.0-1
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
* Tue Sep 23 2008 Julien Dutheil <jdutheil@birc.au.dk> 0.3.0-1
- Initial spec file.
