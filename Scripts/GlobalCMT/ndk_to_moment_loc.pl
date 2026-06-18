#!/usr/bin/env perl
use strict;
use warnings;

sub is_leap_year {
    my ($year) = @_;
    return 1 if $year % 400 == 0;
    return 0 if $year % 100 == 0;
    return $year % 4 == 0;
}

sub day_of_year {
    my ($year, $month, $day) = @_;
    my @month_days = (0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);
    $month_days[2] = 29 if is_leap_year($year);

    my $jday = $day;
    for my $idx (1 .. $month - 1) {
        $jday += $month_days[$idx];
    }
    return $jday;
}

my $line_number = 0;
my ($year, $month, $day, $hour, $minute, $second, $latitude, $longitude, $depth);
my $exponent;

while (my $line = <>) {
    chomp $line;
    my $block_line = $line_number % 5;

    if ($block_line == 0) {
        $line =~ s/[\/:]/ /g;
        my @fields = split ' ', $line;
        ($year, $month, $day, $hour, $minute, $second, $latitude, $longitude, $depth) =
            @fields[1 .. 9];
    }
    elsif ($block_line == 3) {
        my @fields = split ' ', $line;
        $exponent = $fields[0];
    }
    elsif ($block_line == 4) {
        my @fields = split ' ', $line;
        my $moment = $fields[10] * 10 ** $exponent;
        my $jday = day_of_year($year, $month, $day);
        printf "%d %03d %02d %02d %05.2f %g %g %g %g\n",
            $year - 1900, $jday, $hour, $minute, $second,
            $longitude, $latitude, $depth, $moment;
    }

    $line_number++;
}
