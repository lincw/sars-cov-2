#!/usr/bin/env -perl

# use strict;
use warnings;
# use Data::Compare;

open my $file, "<:encoding(utf8)", @ARGV;
my $dummy = <$file>;

my (%mers, %sars1, %sars2, %hku1, %oc43, %h229e, %nl63, @data, $data);
my $count = 0;

while (<$file>) {
    chomp($_);
    @data = split(/,/, $_);
    if ($data[0] eq "MERS") {
        !$mers{$data[2]} ? @{$mers{$data[2]}} = $data[1]: push @{$mers{$data[2]}}, $data[1];
    } elsif ($data[0] eq "HKU1") {
        !$hku1{$data[2]} ? @{$hku1{$data[2]}} = $data[1]: push @{$hku1{$data[2]}}, $data[1];
    } elsif ($data[0] eq "OC43") {
        !$oc43{$data[2]} ? @{$oc43{$data[2]}} = $data[1]: push @{$oc43{$data[2]}}, $data[1];
    } elsif ($data[0] eq "229E") {
        !$h229e{$data[2]} ? @{$h229e{$data[2]}} = $data[1]: push @{$h229e{$data[2]}}, $data[1];
    } elsif ($data[0] eq "SARS1") {
        !$sars1{$data[2]} ? @{$sars1{$data[2]}} = $data[1]: push @{$sars1{$data[2]}}, $data[1];
    } elsif ($data[0] eq "SARS2") {
        !$sars2{$data[2]} ? @{$sars2{$data[2]}} = $data[1]: push @{$sars2{$data[2]}}, $data[1];
    } else {
        !$nl63{$data[2]} ? $nl63{$data[2]} = $data[1]: push @{$nl63{$data[2]}}, $data[1];
    }
    # last if ++$count == 5;
}
close $file;

my @hashes_name = ("MERS", "SARS1", "SARS2", "OC43", "H229E", "NL63", "HKU1");

# for ($i = 0; $i < 7; $i++) {
#     for ($j = 0; $j < 7; $j++) {
my $i = 5;
my $j = 6;
my %old_hash = %nl63;
my %new_hash = %hku1;

        print "$hashes_name[$i] vs $hashes_name[$j]\n";
        print "virus protein\tunique target (sum of 2 virus)\tidentical pair\tdifferent pair\n";
        foreach my $keys (keys %old_hash) {
            my (@union, @intersection, @difference, $union_count);
            my %count = ();

            if ($new_hash{$keys}) {
                foreach my $value ( @{ $old_hash{$keys} } ) {
                    !$count{$value} ? $count{$value} = 1 : $count{$value} ++;
                    $union_count++;
                    }
                foreach my $value ( @{ $new_hash{$keys} } ) {
                    $count{$value} ++;
                    $union_count++;
                }
                foreach $element (keys %count) {
                    push @union, $element;
                    push @{ $count{$element} > 1 ? \@intersection : \@difference }, $element;
                }
                print "$keys\t", $union_count, "\t", scalar @intersection, "\t", scalar @difference, "\n";
            } else {
                next;
            }
        }
#     }
# }
