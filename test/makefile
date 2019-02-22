BINDIR=bin

all: $(BINDIR)/rd_tests

$(BINDIR)/rd_tests: src/main.cpp | $(BINDIR)
	$(CXX) -o $(BINDIR)/rd_tests -Ilib/catch2/single_include/ src/main.cpp

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -rf $(BINDIR)
