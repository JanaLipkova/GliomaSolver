#
#    PSthreads and UthLib thread libraries
#    == Copyright (C) 2007 the OMPi/PSthreads Team
#    == Department of Computer Science, University of Ioannina
#
#    This file is part of PSthreads.
#
#    PSthreads is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    PSthreads is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with PSthreads; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
##   This file must be processed by automake

# Create the default directories if needed

scriptbindir=$(DESTDIR)$(bindir)

install-exec-hook:
	if [ ! -d $(scriptbindir) ]; then \
		mkdir -p $(scriptbindir) ; \
	fi ; \
	cp -f torc_cflags $(scriptbindir); \
	cp -f torc_libs $(scriptbindir);

# Uninstall default directories
uninstall-hook:
	rm -rf $(scriptbindir)/torc_cflags; \
	rm -rf $(scriptbindir)/torc_libs; \
	rm -rf $(scriptbindir);

