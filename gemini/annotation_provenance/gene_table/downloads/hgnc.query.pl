#perl script for downloading specific columns from HGNC
# alternate url for http://tinyurl.com/lvzsywl 

#!/usr/bin/perl -w
use strict;
use LWP::Simple;
my $url = "http://www.genenames.org/cgi-bin/hgnc_downloads?".
    "col=gd_hgnc_id&col=gd_app_sym&col=gd_prev_sym&col=gd_aliases&col=gd_pub_ensembl_id&".
    "status=Approved&status_opt=2&where=&order_by=gd_pub_chrom_map_sort&format=text&".
    "limit=&submit=submit";
my $page = get($url);
print $page;
