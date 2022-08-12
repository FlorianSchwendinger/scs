
docs:
	R -e 'library("roxygen2"); roxygenize(".", c("rd", "namespace"), load_source)'

build:
	R CMD build .

inst: build
	R CMD INSTALL scs*.tar.gz
	
check: build
	R CMD check scs*.tar.gz

cran_check: build
	R CMD check --as-cran scs*.tar.gz

manual: clean
	R CMD Rd2pdf --output=Manual.pdf .

clean:
	rm -f Manual.pdf README.knit.md README.html
	rm -rf .Rd2pdf*

check_mac_m1: build
	R -e "rhub::check(dir(pattern = 'scs_.*.tar.gz'), platform = 'macos-m1-bigsur-release')"

check_mac_old: build
	R -e "rhub::check(dir(pattern = 'scs_.*.tar.gz'), platform = 'macos-highsierra-release-cran')"

check_gcc_san: build
	R -e "rhub::check(dir(pattern = 'scs_.*.tar.gz'), platform = 'linux-x86_64-rocker-gcc-san')"

check_debian_clang: build
	R -e "rhub::check(dir(pattern = 'scs_.*.tar.gz'), platform = 'debian-clang-devel')"

