srcdir = .
prefix = /usr/local
bindir = $(prefix)/bin
mandir = $(prefix)/share/man
datadir = $(prefix)/share
plugindir = $(datadir)/games/plugins
fairydir = $(datadir)/games/fairymax

CC?=gcc
CFLAGS?= -O2 -s
INI_Q?=$(fairydir)/qmax.ini
VERSION?=`grep 'define VERSION' fairymax.c | sed -e 's/.*"\(.*\)".*/\1/'`

ALL= fairymax shamax maxqi fairymax.6.gz

all: ${ALL}

fairymax: fairymax.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -DFAIRYDIR=\"${fairydir}\" fairymax.c -o fairymax

shamax: fairymax.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -DFAIRYDIR=\"${fairydir}\" -DSHATRANJ fairymax.c -o shamax

maxqi: maxqi.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -DINI_FILE=\"${INI_Q}\" maxqi.c -o maxqi

install: ${ALL} ${srcdir}/data/*
	install -d -m0755 $(DESTDIR)$(bindir)
	cp ${srcdir}/fairymax ${srcdir}/shamax ${srcdir}/maxqi $(DESTDIR)$(bindir)
	install -d -m0755 $(DESTDIR)$(fairydir)
	cp ${srcdir}/data/*.ini $(DESTDIR)$(fairydir)
	cp ${srcdir}/data/*.hash $(DESTDIR)$(fairydir)
	install -d -m0755 $(DESTDIR)$(mandir)/man6
	cp ${srcdir}/fairymax.6.gz $(DESTDIR)$(mandir)/man6
	install -d -m0755 $(DESTDIR)$(plugindir)/logos
	cp ${srcdir}/data/logo.png $(DESTDIR)$(plugindir)/logos/fairymax.png
	cp ${srcdir}/data/logo.png $(DESTDIR)$(plugindir)/logos/shamax.png
	cp ${srcdir}/data/logo.png $(DESTDIR)$(plugindir)/logos/maxqi.png
	install -d -m0755 $(DESTDIR)$(plugindir)/xboard
	cp ${srcdir}/data/*.eng $(DESTDIR)$(plugindir)/xboard

fairymax.6.gz: fairymax.pod
	pod2man -s 6 fairymax.pod | gzip -9n > fairymax.6.gz

clean:
	rm -f ${ALL}

dist-clean:
	rm -f ${ALL} *~ data/*~ *.man md5sums

dist: fairymax
	install -d -m0755 Fairy-Max
	install -d -m0755 Fairy-Max/data
	rm -f fairymax.tar fairymax.tar.gz
	cp fairymax.c maxqi.c fairymax.pod Makefile README changelog copyright Fairy-Max
	cp data/* Fairy-Max/data
	(md5sum Fairy-Max/* Fairy-Max/data/* > Fairy-Max/md5sums) || true
	tar -cvvf fairymax-$(VERSION).tar Fairy-Max
	gzip fairymax-$(VERSION).tar
	rm fairymax
	rm Fairy-Max/data/*
	rmdir Fairy-Max/data
	rm Fairy-Max/*
	rmdir Fairy-Max

uninstall:
	rm -f $(DESTDIR)$(plugindir)/logos/fairymax.png
	rm -f $(DESTDIR)$(plugindir)/logos/shamax.png
	rm -f $(DESTDIR)$(plugindir)/logos/maxqi.png
	rm -f $(DESTDIR)$(plugindir)/xboard/fairymax.eng
	rm -f $(DESTDIR)$(plugindir)/xboard/shamax.eng
	rm -f $(DESTDIR)$(plugindir)/xboard/maxqi.eng
	rm -f $(DESTDIR)$(fairydir)/*
	rmdir $(DESTDIR)$(fairydir)
	rm -f $(DESTDIR)$(mandir)/man6/fairymax.6.gz
	rm -f $(DESTDIR)$(bindir)/fairymax
	rm -f $(DESTDIR)$(bindir)/shamax
	rm -f $(DESTDIR)$(bindir)/maxqi

