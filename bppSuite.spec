%define name bppsuite
%define version 0.6.2
%define release 1
%define _prefix /usr

Summary: The Bio++ Program Suite.
Name: %{name}
Version: %{version}
Release: %{release}
Vendor: Julien Dutheil
Source: http://download.gna.org/bppsuite/%{name}-%{version}.tar.gz
License: CeCILL 2
Group: System Environment/Libraries
BuildRoot: %{_builddir}/%{name}-root
Packager: Julien Dutheil
Prefix: %{_prefix}
AutoReq: yes
AutoProv: yes

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

%prep
%setup -q

%build
CFLAGS="-I%{_prefix}/include $RPM_OPT_FLAGS"
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=%{_prefix}"
if [ %{_lib} == 'lib64' ] ; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DLIB_SUFFIX=64"
fi
cmake $CMAKE_FLAGS .
make
make info

%install
rm -rf $RPM_BUILD_ROOT
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
%{_prefix}/share/info/bppsuite.info.gz
%{_prefix}/share/man/man1/bppml.1.gz
%{_prefix}/share/man/man1/bppseqgen.1.gz
%{_prefix}/share/man/man1/bppancestor.1.gz
%{_prefix}/share/man/man1/bpppars.1.gz
%{_prefix}/share/man/man1/bppdist.1.gz
%{_prefix}/share/man/man1/bppconsense.1.gz
%{_prefix}/share/man/man1/bppseqman.1.gz
%{_prefix}/share/man/man1/bppreroot.1.gz
%{_prefix}/share/man/man1/bppphysamp.1.gz
%{_prefix}/share/man/man1/bpptreedraw.1.gz

%changelog
* Thu Jun 09 2011 Julien Dutheil <julien.dutheil@univ-montp2.fr>
- BppSuite 0.6.2 release
* Mon Feb 28 2011 Julien Dutheil <julien.dutheil@univ-montp2.fr>
- BppSuite 0.6.1 release
* Mon Feb 07 2011 Julien Dutheil <julien.dutheil@univ-montp2.fr>
- BppSuite 0.6.0 release
* Thu Mar 25 2010 Julien Dutheil <julien.dutheil@univ-montp2.fr>
- BppSuite 0.5.0 release
* Wed Jun 10 2009 Julien Dutheil <jdutheil@birc.au.dk>
- BppSuite 0.4.0 release
* Thu Dec 11 2008 Julien Dutheil <jdutheil@birc.au.dk>
- BppSuite 0.3.1 release
* Thu Sep 23 2008 Julien Dutheil <jdutheil@birc.au.dk>
- BppSuite 0.3.0 release

