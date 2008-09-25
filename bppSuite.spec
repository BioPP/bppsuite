%define name bppsuite
%define version 0.3.0
%define release 1
%define _prefix /usr/local

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
Requires: Bpp-Utils = 1.3.0
Requires: Bpp-NumCalc = 1.5.0
Requires: Bpp-Seq = 1.4.1
Requires: Bpp-Phyl = 1.6.0

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

%prep
%setup -q

%build
CFLAGS="-I%{_prefix}/include $RPM_OPT_FLAGS" ./configure --prefix=%{_prefix}
make

%install
rm -rf $RPM_BUILD_ROOT
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig

%files
%defattr(-,root,root)
%doc AUTHORS COPYING INSTALL NEWS README ChangeLog
%{_prefix}/bin/bppml
%{_prefix}/bin/bppseqgen
%{_prefix}/bin/bppancestor
%{_prefix}/bin/bppdist
%{_prefix}/bin/bpppars
%{_prefix}/bin/bppseqman
%{_prefix}/bin/bppconsense
%{_prefix}/bin/bppphysamp
%{_prefix}/bin/bppreroot

%changelog
* Thu Sep 23 2008 Julien Dutheil <jdutheil@daimi.au.dk>
- BppSuite 0.3.0 release

