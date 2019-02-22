BINDIR=bin

all: $(BINDIR)/rd

$(BINDIR)/rd: src/main.cpp | $(BINDIR)
	$(CXX) src/main.cpp -o $(BINDIR)/rd

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -rf $(BINDIR)
