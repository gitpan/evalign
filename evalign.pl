#!/usr/bin/perl
# Last Update by /Perl/update_subroutines.pl: Wed Feb 26 15:09:25 GMT 1997
#___________________________________________________________________________
# Title    : evalign.pl
# Function : When you align any sequences by computer algorithms, you want
#            to know whether they are correctly aligned in terms of structures.
#            If the sequences are from already known structures, you can compare
#            and align structural sequences which can be said 'biologically correct'.
#            This program, 'evalign.pl' is for comparing the two sets of sequences
#            aligned, by calculating the absolute position differences between the
#            correct and computer aligned one. This is aware of gap intertions and
#            correct alignment made after wrong alignment segment is counted as correct.
#            It accepts two sequence files at prompt to calculate the differences
#            of positions of the sequences in the input files. The input sequences
#            should be identical in both files.
#            As an option, this also displays Percentage IDentity.
# Usage    : "evalign.pl any_seq_file.msf any_struc_file.jp ["  while any_seq_file.msf
#            is a computer aligned output and any_struc_file.jp is a any seq file
#            from known structures. (eg,  evalign.pl  aa.msf aa.jp )
# Example  : evalign.pl aa.msf aa.jp -ss -H -E -p
#
# Argument : Two files of sequence alignment. The first one should be COMPUTER aligned
#            and the second one is the CORRECT (i.e., structural) alignment.
# Options  : seg is for showing the accuracy of alignment on secondary str. blocks.
#            ss  is for showing DSSP secondary structure assignment in output.
#            H   is for showing HELIX DSSP secondary structure assignment in output.
#            E   is for showing Beta-strand DSSP secondary structure assignment in out.
#            s   is for sorted final output.
#            p   is for displaying conventional percent ID.
#            h   is for displaying help
#            ns  is for $no_simplify by -ns, ns, Ns, NS, -Ns # seq names are sorted in final output
#            t=  is for convert to num of 1 or 0 threshold.
#            c   is for convert to num of 1 or 0, default threshold '1' is used
#            N   for DO NOT Normalize the error rate which can be more than 1 digit
#
# $NO_normalize      = 1  by  N -N
# $segment_rate      = 1  by  -seg, seg, Seg # Shows secondary str. block PSR
# $show_percent_id   = 1  by  -p, -P, p, P,  # Shows conventional percent ID.
# $show_sec_str      = 1  by  -ss, ss or SS  # Show Secondary Structure -ss option
# $HELIX_only        = H  by  -H, H          # Shows conventional percent ID.
# $BETA_only         = E  by  -E, E          # Shows conventional percent ID.
# $print_sort        = s  by  -s, s or S     # seq names are sorted in final output
# $interlaced        = i  by  -i, i or I     # seq names are sorted in final output
# $no_simplify       = 1  by  -ns, ns, Ns, NS, -Ns # seq names are sorted in final output
# $threshold         =    by  t=    # seq names are sorted in final output
# $convert_to_0_or_1 = 1  by  -c, C, c, Con  # seq names are sorted in final output
# $HELP              = 1  by  -h, h          # for showing help
#
# Returns  : simple shifted positions.
# Author   : A Biomatic
# Version  : 1.4
# Package  : Part of Bioperl project.
# Warning  :
#-------------------------------------------------------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
# Start of the main program
#___________________________________
   my ($file1, $file2)=@{&parse_arguments(2)};
   my($real_file_name, $n, $m, @dssp_hash, @names_of_seqs,%result, %sec_hash1,
	  %big_hash_sec_str, %sec_hash2 );
   my($group_name, $extension)=split(/\./, $file1);
   my(%array1) =%{&read_any_seq_files(\$file1)};
   my(%array2) =%{&read_any_seq_files(\$file2)};

   &show_hash(\%array1) if $verbose==1 or $verbose eq 'v';

   @names_of_seqs = sort keys %array1;
   #  @names_of_seqs has 1cdg, 6taa, 2aaa, etc.

   print "\n# \@names_of_seqs are: @names_of_seqs\n"  if $verbose==1 or $verbose eq 'v';

   %result  = %{&get_position_shift_rate(\%array1, \%array2)};
   %result  = %{&normalize_numbers(\%result)} unless $NO_normalize==1;

   $convert_to_0_or_1 = 1 if defined($threshold);
   $threshold=1 unless defined($threshold);
   if($convert_to_0_or_1 ==1){
	 %result  = %{&convert_num_0_or_1_hash_opposite(\%result, \$threshold )};
   }

   ##############################################################
   ###   When "Show Secondary Structure" '-ss' option is used ###
   ##############################################################
   if( $show_sec_str == 1){
	   for $name( @names_of_seqs ){
		  my(%temp_single_hash);
		  if( defined(open_dssp_files2) ){
			 (@out_ref) = @{&open_dssp_files2(\$name, \$HELIX_only, \$BETA_only)}; }
		  else{ (@out_ref) = @{&open_dssp_files(\$name, \$HELIX_only, \$BETA_only)}; }

		  show_hash(@out_ref); ## Checking if open_dssp_files succeeded

		  %{"$name"} =%{$out_ref[0]};  ### %{"$name"} is one hash of one key and one Value.
									   ### like '1cdg HHHHHHHHHHHHHH   EEEEEEEE EE HHHHH'.
		  $temp_single_hash{$name}=$array2{$name};  ##  $array2{$name} has the structural sequence of
													## one seq name.
		  %{"$name"}=%{&insert_gaps_in_seq_hash(\%{"$name"}, \%temp_single_hash)};
		  push(@dssp_hash, \%{"$name"});
	   }   ## this is the returning array of ref.
	   #""""""""""""""""""""""""""""""""""""""""""""""""""""""
	   ### If show  wrong segment rate  option -seg is set ###
	   #""""""""""""""""""""""""""""""""""""""""""""""""""""""
	   if($segment_rate == 1){
		  for ($n=0; $n < @dssp_hash; $n ++ ){
			 my(%combined);
			 %sec_hash1 = %{$dssp_hash[$n]};
			 for($m=($n+1); $m < @dssp_hash; $m ++ ){
				%sec_hash2 = %{$dssp_hash[$m]};
				%combined = %{&overlay_seq_for_identical_chars(\%sec_hash1, \%sec_hash2, '-')};
				%big_hash_sec_str = (%combined, %big_hash_sec_str); ## adding up single hashes
																	## to bigger one.
			 }
		  }
	   }
   }
   if( $show_seq_align==1){
	  &print_seq_in_block(\%result, \%array2, @dssp_hash, @combined_dssp_hash,
						  \$interlaced, \$print_sort,  50);
   }else{
	  &print_seq_in_block(\%result, @dssp_hash, @combined_dssp_hash,
						  \$interlaced, \$print_sort,  50);
   }

   #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
   #      When "Show percent identity" '-p' option is used
   #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
   if($show_percent_id ==1){
	 my(@identities) = @{&get_correct_percent_alignment_rate(\$file1, \$file2)};
	 print "\n", "_"x88, "\n";
	 print "\nFinal output:  @identities\n";
	 print "\n", "_"x88, "\n";
   }





#____________________________________________________________________
# Title    : get_segment_shift_rate
# Function : calculates the secondary structure segment shift rate.
# Usage    : &get_segment_shift_rate(\%hash_for_errors, \%hash_for_sec_str);
# Example  : <input example> First block is for the first hash input
#                            and Second is for the second hash input.
#
#           1cdg_6taa      00000442222222222242222222222777700000007000000000
#           1cdg_2aaa      00000442222222222242222222222777700000007000000000
#           2aaa_6taa      00000000000000000000000000000000000000000000000000
#
#           1cdg_6taa      -------EEE-----------EE--EEEE------EE---------EEE-
#           1cdg_2aaa      -------EEE-----------EE--EEEE------EE---------EEE-
#           2aaa_6taa      -------EEEEE------EE-EEEEEEEE----EEEE-------EEEEE-
#
#            <intermediate output example>
#           2aaa_6taa      -------00000---------00000000----0000-------00000-
#           1cdg_6taa      -------442---------------2222-----------------000-
#           1cdg_2aaa      -------222---------------2222-----------------000-
#
#            <Final output>
#           2aaa_6taa      0%
#           1cdg_6taa      67%
#           1cdg_2aaa      67%
#
# Argument :
# Returns  :
# Options  : 'p' or 'P' for percentage term(default)
# Options  : 'r' or 'R' for ratio term (0.0 - 1.0), where 1 means all the
#             segments were wrongly aligned.
# Options  : 's' or 'S' for Shift rate (it actually caculates the position shift
#             rate for the secondary structure segment.
# Options  : 'h' or 'H' for position Shift rate (it actually caculates the position
#             shift rate for helical segments). If this is the only option, it
#             will show the default percentage term rate for helical segments.
#             If used with 'r', it will give you ratio (0.0 - 1.0) for helical
#             segment. If used with 's' option, it will give you position shift
#             rate for only helical segments.
# Options  : 'e' or 'E' for position Shift rate (it actually caculates the position
#             shift rate for beta segments). If this is the only option, it will
#             show the default percentage term rate for beta segments. If used
#             with 'r', it will give you ratio (0.0 - 1.0) for beta. If used
#             with 's' option, it will give you position shift rate for only
#             beta segments.
# Tips     :
# Author   : A Biomatic
# Version  : 1.0
# Package  :
# Keywords : later sub of get_position_shift_rate for secondary structure regions
#            get positon shift rate for secondary structure regions.
# Warning  :
#---------------------------------------------------------------
sub get_segment_shift_rate{
  my($i, $k, $j, @hash_input, $option_string, %h, %superposed_hash,
	 $name, %out, $gap_chr, @str1, @str2, %temp, %hash_error, %hash_secondary);
  print chr(7) , "\n chr(7) get_segment_shift_rate is called \n";
  ##########################################
  #       general argument handling        #
  ##########################################
  for($k=0; $k< @_ ;$k++){
	 if( ( !ref($_[$k]) )&&($_[$k]=~ /^(\w)$/) ){
		$option_string  .= $1;    }
	 elsif((ref($_[$k]) eq "SCALAR")&&(${$_[$k]}=~ /^(\w)$/) ){
		$option_string  .= $1;    }
	 elsif(ref($_[$k]) eq "HASH") {
		%temp = %{$_[$k]};
		my(@keys)= sort keys (%temp);
		my($temp_seq) = $temp{$keys[0]};

		if($temp_seq=~/\d\d+/){
		   %hash_error = %temp; }
		else{ %hash_secondary = %temp; }
	 }
  }# """""""" OUTPUT are  : %hash_error  &  %hash_secondary
  ##########################################
  #       general argument handling end    #
  ##########################################

  %hash_secondary =%{&tidy_secondary_structure_segments(\%hash_secondary)};
  print_seq_in_block(\%hash_error);
  print_seq_in_block(\%hash_secondary);
  %superposed_hash =%{&superpose_seq_hash(\%hash_error, \%hash_secondary)};
  print_seq_in_block(\%superposed_hash);
  %count_out = %{count_sequence_segments(\%superposed_hash)};
  %h=%{&get_wrong_segment_rate(\%superposed_hash)};
  \%h;
}


#________________________________________________________________________
# Title     : overlay_seq_by_certain_chars
# Usage     : %out =%{&overlay_seq_by_certain_chars(\%hash1, \%hash2, 'HE')};
# Function  : (name1 000000112324)+(name1  ABC..AD..EFDK ) => (name1 000..00..12324)
#             (name2 000000112324)+(name2  --HHH--EEEE-- ) => (name1 ---000--1123--)
#             uses the second hash a template for the first sequences. gap_char is
#             '-' or '.' or any given char or symbol.
#             To insert gaps rather than overlap, use insert_gaps_in_seq_hash
# Example   : %out =%{&overlay_seq_by_certain_chars(\%hash1, \%hash2, 'E')};
#             output> with 'E' option >>> "name1     --HHH--1232-"
# Warning   : If gap_chr ('H',,,) is not given, it replaces all the
#             non-gap chars (normal alphabet), ie,
#             it becomes 'superpose_seq_hash'
# Class     :
# Keywords  : Overlap, superpose hash, overlay, superpose_seq_hash
# Options   : E for replacing All 'E' occurrances in ---EEEE--HHHH----, etc.
#             : H for replacing all 'H'  "     " "
# Package   : Array_Util
# Reference :
# Returns   : one hash ref.
# Tips      :
# Argument  : 2 ref for hash of identical keys and value length.
# Todo      :
# Author    : A Biomatic
# Version   : 1.0
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub overlay_seq_by_certain_chars{
  my($i, $k,$j, $name, @in, %out, $gap_chr, @str1, @str2);
  ######################################
  ####### Sub argument handling ########  $gap_chr here can be 'HE' etc.
  ######################################
  for($k=0; $k< @_ ;$k++){
	  if( ( !ref($_[$k]) )&&($_[$k]=~ /^(.+)$/) ){
		  $gap_chr  .= $1;
	  }elsif((ref($_[$k]) eq "SCALAR")&&(${$_[$k]}=~ /^(.)+$/) ){
		  $gap_chr  .= $1;
	  }elsif(ref($_[$k]) eq "HASH") { push(@in,  $_[$k]); }
  }

  if($#in < 1){
	  print "\n overlay_seq_by_certain_chars needs 2 hashes. Error \n"; exit; }
  my(%hash1)=%{$in[0]};
  my(%hash2)=%{$in[1]};
  my(@names1)= sort keys %hash1;
  my(@names2)= sort keys %hash2;
  (@names1 > @names2)? $bigger=@names1 : $bigger=@names2;
  for ($j=0; $j < $bigger; $j++){
	  @str1=split(//, $hash1{$names1[$j]});
	  @str2=split(//, $hash2{$names2[$j]});
	  if( ($gap_chr eq '') && ($hash2{$names2[$j]}=~/(\W)/) ){
		  $gap_chr=$1;
		  for($i=0; $i < @str2; $i++){
			  if($str2[$i] =~ /$gap_chr/){ $str1[$i]=$gap_chr;}     }
		  $out{$names1[$j]}=join(",", @str1);
	  }else{
		  for($i=0; $i < @str2; $i++){
			  if($gap_chr =~ /$str2[$i]/){ $str2[$i]=$str1[$i];}    }
		  $out{$names1[$j]}=join(",", @str2);    }
  }
  \%out;
}




#________________________________________________________________________
# Title     : read_dir_names_only
# Usage     : @all_dirs_list = @{&read_dir_names_only(\$absolute_path_dir_name, ....)};
# Function  : read any dir names and and then put in array.
# Example   :
# Warning   : This does not report '.', '..'
#             Only file names are reported. Compare with &read_any_dir
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : one ref. of array.
# Tips      :
# Argument  : takes one or more scaler references. ('.', \$path, $path, ... )
# Todo      :
# Author    : A Biomatic
# Version   : 3.1
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub read_dir_names_only{
  my($in_dir, $i,$k, @possible_dirs,
	  @final_files, $full_dir, $pwd, $path,@read_files);
  $pwd=`pwd`; chomp($pwd); $full_dir=1;
  for($k=0; $k < @_; $k++){
	 if   ( ($_[$k] eq '.') || !(defined($_[$k]))){  $in_dir=$pwd;  }
	 elsif(!(ref($_[$k]))){   $in_dir=$_[$k];   }
	 elsif(ref($_[$k])){      $in_dir =${$_[$k]};    }
	 if($in_dir =~ /^([\w\-\.]+)$/){  $in_dir="$pwd\/$in_dir"; $full_dir = 0; }
	 else{ $full_dir =1; }
	 ##########  Main READING PART ##########
	 opendir(DIR1,"$in_dir");
	 @read_files = readdir(DIR1);
	 for($i=0; $i < @read_files; $i ++){
		$read_files[$i]="$in_dir\/$read_files[$i]";
		if( ($read_files[$i] !~ /\/\.\.?$/) && ( -d $read_files[$i]) ){
		  $read_files[$i]=~s/\.\///; ## removing ./ in front of dirs (in bash)
		  push(@final_files, "$read_files[$i]");
		}
	 }
  }
  return([sort @final_files]);
}

#________________________________________________________________________
# Title     : overlay_seq_for_identical_chars
# Usage     : %out =%{&overlay_seq_for_identical_chars(\%hash1, \%hash2, '-')};
# Function  : (name1         --EHH--HHEE-- )
#             (name2         --HHH--EEEE-- ) ==> result is;
#
#             (name1_name2   -- HH--  EE-- )
#             to get the identical chars in hash strings of sequences.
#
# Example   : %out =%{&overlay_seq_for_identical_chars(\%hash1, \%hash2, '-')};
#             output> with 'E' option >>> "name1     --HHH--1232-"
# Warning   : Works only for 2 sequence hashes.
# Class     :
# Keywords  : Overlap, superpose hash, overlay identical chars, superpose_seq_hash
# Options   :
# Package   : Array_Util
# Reference :
# Returns   : one hash ref. of the combined key name (i.e., name1_name2). Combined by '_'
# Tips      :
# Argument  : 2 ref for hash of identical keys and value length. One optional arg for
#             replacing space char to the given one.
# Todo      :
# Author    : A Biomatic
# Version   : 1.0
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub overlay_seq_for_identical_chars{
  my($i, $k,$j, $name1, $name2, @in, %out, @out_chars, $gap_chr, @str1, @str2);
  ######################################
  ####### Sub argument handling ########  $gap_chr here can be 'HE' etc.
  ######################################
  for($k=0; $k< @_ ;$k++){
	  if( ( !ref($_[$k]) )&&($_[$k]=~ /^(.)$/) ){
		  $gap_chr  .= $1;    }
	  elsif((ref($_[$k]) eq "SCALAR")&&(${$_[$k]}=~ /^(.)$/) ){
		  $gap_chr  .= $1;    }
	  elsif(ref($_[$k]) eq "HASH") { push(@in,  $_[$k]); }    }

  if(@in < 2){ print "\n overlay_seq_for_identical_chars needs 2 hashes. Error \n"; exit; }
  my(%hash1)=%{$in[0]}; my(%hash2)=%{$in[1]};
  my(@names1)=sort keys %hash1; my(@names2)= sort keys %hash2;
  $name1 = $names1[0]; $name2 = $names2[0];
  @str1=split(/|\,/, $hash1{$names1[0]}); @str2=split(/|\,/, $hash2{$names2[0]});
  for($i=0; $i < @str1; $i++){
	  if($str1[$i] eq $str2[$i] ){
		  push(@out_chars, $str1[$i]); }
	  elsif( defined($gap_chr) ){ push(@out_chars, $gap_chr); }
	  else{ push(@out_chars, ' '); }
  }
  if( $name2 gt $name1){
	  $out{"$name1\_$name2"}= join(",", @out_chars); }
  else{  $out{"$name2\_$name1"}= join(",", @out_chars); }
  \%out;
}

#________________________________________________________________________
# Title     : normalize_numbers ( from any numbers to  0 - 9 )
# Usage     : %output=%{&normalize_numbers(\%hash1)};
#             originally made to normalize the result of get_posi_rates_hash_out
#             in   'scan_compos_and_seqid.pl'
# Function  : with given numbers in hashes, it makes a scale of 0-9 and puts
#             all the elements in the scale. Also returns the average of the numbs.
# Example   : intputhash>                   Outputhash>
#             ( '1-2', '12,.,1,2,3,4',     ( '1-2',   '9,.,0,1,2,3',
#              '2-3', '12,.,1,5,3,4',       '2-3',   '9,.,0,4,2,3',
#              '4-3', '12,3,1,2,3,4',       '3-1',   '9,3,.,.,2,3',
#              '3-1', '12,4,.,.,3,4' );     '4-3',   '9,2,0,1,2,3' );
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : (\%norm_hash1, \%norm_hash2, \%norm_hash3,.... )
#
# Tips      :
# Argument  : (\%hash1, %hash2, \%hash3, ....)
# Todo      :
# Author    : A Biomatic
# Version   : 1.0
# Used in   : evalign.pl
# Enclosed  :
#--------------------------------------------------------------------
sub normalize_numbers{
  my(@in)=@_;
  my($split_char)=',';
  my(@out_ref_of_hash, $min, $max, $name, $u,$sum, $av, $range, @num_array,%in);
  ($min, $max, $sum, $av)=&main::hash_stat_for_all(@in);
  if(($max-$min)==0){ $range = 1} else { $range= ($max -$min) };
  for ($u=0; $u< @in ; $u++){
	  %in=%{$_[$u]};
	  my(@keys) = keys %in;
	  if($in{$keys[0]}=~/\,/){ $split_char=','; }  else { $split_char=''; };
	  for $name (@keys){  @num_array = split(/$split_char/, $in{$name});
		  for (@num_array){ $_ = int(($_ / $range)*8) if ($_ =~ /[\-]*\d+/); }
		  $in{$name}=join("$split_char", @num_array); }
	  push(@out_ref_of_hash, \%in);  }
  if( @out_ref_of_hash == 1)  {  return( $out_ref_of_hash[0]); }
  elsif( @out_ref_of_hash > 1){  return( @out_ref_of_hash   ); }
}

#________________________________________________________________________
# Title     : hash_stat_for_all
# Usage     : %out=%{&hash_average(\%in, \%in2,..)};
# Function  : gets the min, max, av, sum for the whole values of ALL the
#             hashes put in. (grand statistics)
# Example   : %in =(1, "13242442", 2, "92479270", 3, "2472937439");
#             %in2=(1, "28472", 2, "23423240", 3, "123412342423439");
#
#             %in =(name1, "1,3,2,4,2,4,4,2", name2, "9,2,4,7,9,2,7,0");
#
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : normal array of ($min, $max, $sum, $av)
#             Example  out:>                 |  min max sum  av
#                            -----------------------------------
#                            of the whole    |   0   9  110   6
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   : 1.0
# Used in   : normalize_numbers
# Enclosed  :
#--------------------------------------------------------------------
sub hash_stat_for_all{
  package hash_stat_for_all;
  my($elem,@out_av_hash, @out_array,$v,@num_arr,$sum,$min, $av,$num_all, $max,$split_char);
  for($v=0; $v<@_; $v++){
	  my(%input)=%{$_[$v]};
	  for $name(keys(%input)){
		 if($input{$name} =~ /\,/){ $split_char=',';
		 }else{
		    $split_char='';
		 }
		 @num_arr=split(/$split_char/, $input{$name});
		 for $elem(@num_arr){
			if($elem =~/[\-]*\d+/){
			    $min=$elem unless(defined($min));
				$min =$elem if $elem < $min; $max =$elem if $elem > $max;
				$sum+=$elem; $num_all++;
			}
		 }
	  }
  }
  if($num_all == 0){ $av=0; $sum=0; $min=0; $max=0; }
  else { $av=$sum/$num_all; }
  push(@out_array, ($min, $max, $sum, $av));
  package main;
  return(@out_array);
}


#________________________________________________________________________
# Title     : tidy_secondary_structure_segments
# Usage     : print_seq_in_block(&tidy_secondary_structure_segments(\%hash, 'e4', 'h4'), 's');
#
# Function  : receives any secondary structure assignment hashes and
#             tidys up them. That is removes very shoft secondary structure
#             regions like( --HH--, -E-, -EE- ) according to the given minimum
#             lengths(threshold) of segments by you.
# Example   : print_seq_in_block(&tidy_secondary_structure_segments(\%hash, 'e4', 'h4'), 's');
#             <makes following into the next block>
#
#             1cdg_2aaa      -------EEE-----------EE--EEEE------EE---------EEE-
#             1cdg_6taa      -------EEE-----------EE--EEEE------EE---------EEE-
#             2aaa_6taa      -------EEEEE------EE-EEEEEEEE----EEEE-------EEEEE-
#
#             <example output>
#
#             1cdg_6taa      -------------------------EEEE---------------------
#             1cdg_2aaa      -------------------------EEEE---------------------
#             2aaa_6taa      -------EEEEE---------EEEEEEEE----EEEE-------EEEEE-
#
# Warning   :
# Class     :
# Keywords  :
# Options   : something like 'H3' or 'E3' for minimum segment length set to 3 positions.
# Package   : Bio::Seq
# Reference :
# Returns   : array of references of hashes.
# Tips      :
# Argument  : hashes and [options]. No options result in default of 'H3', 'E3'
# Todo      :
# Author    : A Biomatic
# Version   : 1.0.0
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub tidy_secondary_structure_segments{
  my($i, $k,$a, $j, $helix_min, $beta_strand_min, %hash, @keys, @hash,
	  $option_string, @hash_out, $string1, $name, %out, $gap_chr, @str1, @str2,
	  @stringout, @string_segH, @string_segE, $countH, $countE
	  );

  #### Default helix and beta strand segment length setting #####
  $helix_min=3;
  $beta_strand_min=3;

  ########################################################################
  #####   general argument handling  for options of segment length  ######
  ########################################################################
  for($k=0; $k< @_ ;$k++){
	  if( ( !ref($_[$k]) )&&($_[$k]=~ /^[Hh](\d+)$/) ){
		  $helix_min  = $1;    }
	  elsif( ( !ref($_[$k]) )&&($_[$k]=~ /^[Ee](\d+)$/) ){
		  $beta_strand_min  = $1;    }
	  elsif((ref($_[$k]) eq "SCALAR")&&(${$_[$k]}=~ /^[Hh](\d+)$/) ){
		  $helix_min  = $1;    }
	  elsif((ref($_[$k]) eq "SCALAR")&&(${$_[$k]}=~ /^[EeBb](\d+)$/) ){
		  $beta_strand_min  = $1;    }
	  elsif(ref($_[$k]) eq "HASH") { push(@hash,  $_[$k]); }    }

  for($i=0; $i < @hash; $i++){
	  my(%hash) = %{$hash[$i]};
	  @keys = sort keys( %hash );
	  for($j=0; $j < @keys; $j++){
		  my(@string_segH, @string_segE, @stringout);
		  $string1=$hash{$keys[$j]};
		  $gap_char = $1 if ($string1=~ /(\W)/);

		  ##### actual cleaning ####
		  my(@string) = split(//, $string1);
		  for($a = 0; $a < @string; $a++){
			 if($string[$a] !~/[HE]/){ ### if the splited element doesn't match 'H' or 'E'

				 ##### If any of the HH or EE counter is over the given minimum($helix_min,,)
				 if((@string_segH >= $helix_min)||( @string_segE >=$beta_strand_min)){
					 push(@stringout, @string_segH, @string_segE, '-');
					 @string_segH=@string_segE=();     }   ## just resetting.
				 else{  ### if the accumulated 'HH' or 'EE' is smaller than the minimum
					 for(0.. (@string_segH + @string_segE) ){
						push(@stringout, '-'); ### replace the short 'EE' etc with '-'
					 }
					 @string_segH=@string_segE=();  ## just resetting.
				 }
			 }
			 elsif($string[$a] =~ /^([Hh])$/){
				 push(@string_segH, $1); }
			 elsif($string[$a] =~ /^([Ee])$/){
				 push(@string_segE, $1); }
		  }
		  $hash{$keys[$j]}=join("", @stringout);
	  }
	  push(@hash_out, \%hash);
  }
  if(@hash_out == 1){ return($hash_out[0]);
  }elsif(  @hash_out > 1 ){ return(@hash_out); }
}
#________________________________________________________________________
# Title     : superpose_seq_hash
# Usage     : %out =%{&superpose_seq_hash(\%hash1, \%hash2)};
# Function  : (name1 000000112324)+(name1  ABC..AD..EFD ) => (name1 000..01..324)
#             uses the second hash a template for the first sequences. gap_char is
#             '-' or '.'
#             To insert gaps rather than overlap, use insert_gaps_in_seq_hash
# Example   :
# Warning   : Accepts only two HASHes and many possible gap_chr. Default gap is '-'
# Class     :
# Keywords  : overlay sequence, overlay alphabet, superpose sequence,
# Options   :
# Package   :
# Reference :
# Returns   : one hash ref.
# Tips      :
# Argument  : 2 refs. for hash of identical keys and value length and gap_chr.
# Todo      :
# Author    : A Biomatic
# Version   : 1.0
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub superpose_seq_hash{
  if($debug eq 1){ print __LINE__, " # superpose_seq_hash : \n"; }
  my($gap_chr)='-';
  my($i, $j, %hash1, %hash2, $name, %out, @str1, @str2);

  if((ref($_[0]) eq 'HASH')&&(ref($_[1]) eq 'HASH')){
	  %hash1=%{$_[0]}; %hash2=%{$_[1]}; }
  else{ print "\n superpose_seq_hash needs hash ref\n"; print chr(007); exit; }

  my(@names1)=sort keys %hash1; my(@names2)=sort keys %hash2;
  (@names1 > @names2)? $bigger=@names1 : $bigger=@names2;

  for ($j=0; $j < $bigger; $j++){
	 if($hash2{$names2[$j]}=~/(\W)/){ $gap_chr = $1; }
		@str1=split(/|\,/, $hash1{$names1[$j]});
		@str2=split(/|\,/, $hash2{$names2[$j]});
		for($i=0; $i < @str2; $i++){
		  if($str2[$i] ne $gap_chr){ $str2[$i]=$str1[$i];  } }
		$out{$names1[$j]}=join(",", @str2);
  }
  \%out;
}

#________________________________________________________________________
# Title     : insert_gaps_in_seq_hash
# Usage     : %out_extended_seq =%{&insert_gaps_in_seq_hash(\%hash1, \%hash2)};
# Function  : superpose two hashes of the same sequence or same seq. length sequences,
#             but unlike 'superpose_seq_hash', this inserts gaps and extend the
#             sequences.
#             (name1_sec  hHHHHHH EEEEEEE) +
#             (name1_seq  .CDEABC..AD..EFD..EKST) => (name1_ext  .hHHHHH..H...EEE..EEEE)
#             In the example, the undefined sec. str. position is replaced as gaps('.')
#             Uses the second hash a template for the first sequences. gap_char is
#             '-' or '.'
# Example   :
# Warning   : coded by A Biomatic
# Class     :
# Keywords  : superposing sequences with gaps
# Options   :
# Package   :
# Reference :
# Returns   : one hash ref.
# Tips      :
# Argument  : 2 ref for hash of identical keys and value length.
# Todo      :
# Author    : A Biomatic
# Version   : 1.1
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub insert_gaps_in_seq_hash{
  my($gap_char)='-';
  my($i, $j, $t, %hash1, %hash2, $bigger, $name, %out, @str1, @str2);
  my($join_char) =',';  ## <<-- This is for the final output joined by this var.

  if((ref($_[0]) eq 'HASH')&&(ref($_[1]) eq 'HASH')){
	  %hash1=%{$_[0]};
	  %hash2=%{$_[1]};
  }
  else{
	  print "\n superpose_seq_hash needs hash ref\n";
	  print chr(007); exit;
  }
  my(@names1)=keys %hash1;
  my(@names2)=keys %hash2;
  (@names1 > @names2)? $bigger=@names1 : $bigger=@names2;

  for ($j=0; $j < $bigger; $j++){
		if( $hash2{$names2[$j]}=~/(\W)/){
			$gap_char = $1; } ## <<-- finding the gap_char
		$hash1{$names1[$j]}=~ s/ /$gap_char/g; ## <<-- replacing space with 'gap_char';

		@str1=split(/|\,/, $hash1{$names1[$j]});
		@str2=split(/|\,/, $hash2{$names2[$j]});
		for($i=0 ; $i < @str2; $i++){
			if($str2[$i] =~ /\w/){
				$str2[$i] = shift @str1;
			}
		}
		$out{$names1[$j]}=join(",", @str2);
  }
  return(\%out);
}

#________________________________________________________________________
# Title     : get_position_shift_rate
# Usage     : %rate_hash = %{&get_position_shift_rate(\%hash_msf, \%hash_jp)};
# Function  : This is to get position specific error rate for line display rather than
#             actual final error rate for the alignment. Takes two file names of seq.
#             Output >>
#             seq1_seq2  1110...222...2222
#             seq2_seq3  1111....10...1111
#             seq1_seq3  1111....0000.0000
#
# Example   : my(%error_rate)=%{&get_position_shift_rate(\%input, \%input2)};
# Warning   : split and join char is ','; (space)
# Class     :
# Keywords  :
# Options   : 'ss' for secondary structure regions(Helix and Beta region only
#                 calculation for error rate). There is specialized sub called
#              get_segment_shift_rate for sec. str. only handling.
#
#    $ss_opt            becomes    ss by  ss, SS, -ss, -SS     #  for secondary structure only
#    $H                 =         'H' by   -H or -h or H       # to retrieve only H segment
#    $S                 becomes   'S' by   -S or  S            # to retrieve only S segment
#    $E                 becomes   'E' by   -E or  E            # to retrieve only E segment
#    $T                 becomes   'T' by   -T or -t or T or t  # to retrieve only T segment
#    $I                 becomes   'I' by   -I or  I            # to retrieve only I segment
#    $G                 becomes   'G' by   -G or -g or G or g  # to retrieve only G segment
#    $B                 becomes   'B' by   -B or -b or B or b  # to retrieve only B segment
#    $HELP              becomes    1  by   -help   # for showing help
#    $simplify          becomes    1  by   -p or P or -P, p
#    $simplify          becomes    1  by   -simplify or simplify, Simplify SIMPLIFY
#    $comm_col          becomes   'C' by   -C or C or common
#    $LIMIT             becomes    L  by   -L, L               # to limit the error rate to 9 .
#
# Package   :
# Reference :
# Returns   : \%final_posi_diffs;
# Tips      :
# Argument  : %{&get_position_shift_rate(\%msfo_file, \%jpo_file)};
#             Whatever the names, it takes one TRUE structral and one ALIGNED hash.
# Todo      :
# Author    : A Biomatic
# Version   : 1.5
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub get_position_shift_rate{
	#"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
	my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
	my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
	my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
	if($debug==1){print "\n\t\@hash=\"@hash\"
	\@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
	\@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	print __LINE__," \$char_opt is  \"$char_opt\" in get_position_shift_rate\n" if $debug eq 1;
	print __LINE__," \@string is  \"@string\" in get_position_shift_rate\n" if $debug eq 1;
	print __LINE__," \$LIMIT is  \"$LIMIT\" in get_position_shift_rate\n" if $debug eq 1;

	my(%arraySEQ)=%{$hash[0]};
	my(%arraySTR)=%{$hash[1]};
	my($gap_char, %final_posi_diffs, @stringSTR,@stringSEQ,@seq_positionSEQ,
		@seq_positionSTR,$len_of_seq, @position_diffs, @position_corrected1,
		@names, @whole_length, %array3, @keys_common, %DSSP_common, @stringDSSP_common);

	$gap_char='.';

	%arraySTR = %{&hash_common_by_keys(\%arraySTR, \%arraySEQ)};
	%arraySEQ = %{&hash_common_by_keys(\%arraySEQ, \%arraySTR)};
	%arraySEQ = %{&remov_com_column(\%arraySEQ)};
	%arraySTR = %{&remov_com_column(\%arraySTR)};

	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	if($debug eq 1){
		print __LINE__,
		" ## sorting sequence names. To make things constant. \n\n";  }
	@names= sort keys %arraySTR;
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	#  If common column of secondary structure representation option $comm_col is set
	#  open_dssp_files sub routine will get the common seq parts of all the sequences.
	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	if($comm_col =~ /C/i){
		%DSSP_common=%{&open_dssp_files( @names, $H, $S, $E, $T, $I, $G, $B, $simplify, 'C')};
		@keys_common= keys %DSSP_common;
		@stringDSSP_common = split(/|\,/, $DSSP_common{$keys_common[0]});
		if($debug2 eq 1){ print __LINE__," \$comm_col is set to: $comm_col \n";
			print __LINE__," \@stringDSSP_common is :@stringDSSP_common \n";
		}
	}

	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	# Comparing two hashes
	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	for $name (@names){
		#"""""""""""""""" Splitting the sequence string
		if($arraySEQ{$name}=~/\,\S+\,/){
			@stringSEQ =split(/\,/, $arraySEQ{$name});
			@stringSTR=split(/\,/, $arraySTR{$name});  }
		else{
			@stringSEQ =split(//, $arraySEQ{$name});
			@stringSTR=split(//, $arraySTR{$name});
		}
		print "\n",__LINE__, " \@stringSEQ  is  @stringSEQ \n" if $debug2 eq 1;
		print "\n",__LINE__, " \@stringSTR  is  @stringSTR \n" if $debug2 eq 1;

		#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		#   Contracting  the SEQ.
		#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		@seq_positionSEQ = @{&get_posi_sans_gaps(\$arraySEQ{$name})};
		@seq_positionSTR = @{&get_posi_sans_gaps(\$arraySTR{$name})};

		#"""""""""""""""" To get secondary structure only calc  """"""""""""""""""""""""""""
		# It superposes the NON sec. region on  @seq_positionSTR to nullify positions.
		#  get_posi_diff ignores non char positions in calc.
		#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		if( ($ss_opt =~ /ss$/i) && ($comm_col !~ /C/i) ){
			%DSSP=%{&open_dssp_files($name, $H, $S, $E, $T, $I, $G, $B, $simplify, $comm_col)};
			if($debug1 eq 1){
			   print "\n",__LINE__," open_dssp_files has options \$H ->$H \$S->$S \$E->$E \n";
			   print "\n",__LINE__," \$T->$T \$I->$I \$G->$B \$simplify->$simplify \$comm_col ->$comm_col\n";
			   &show_hash( \%DSSP );
			}
			if(ref(\%DSSP) eq 'HASH'){ # to check it %DSSP was valid, If not it skips overlaying
				@stringDSSP = split(/|\,/, $DSSP{$name});
				$size_of_stringDSSP = @stringDSSP;
				$size_of_seq_positionSTR = @seq_positionSTR;
				if($debug2 eq 1){
					  print "\n",__LINE__," \@stringDSSP is \n @stringDSSP\n";
					  print "\n",__LINE__," Size of \@stringDSSP      is $size_of_stringDSSP\n" ;
					  print "\n",__LINE__," Size of \@seq_positionSTR is $size_of_seq_positionSTR\n";
					  print "\n",__LINE__," \$gap_char is \"$gap_char\" \n" ;
				}
				#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
				#   When the sec. str is not defined in DSSP, I delete the position of
				#   @stringDSSP to gap(ie. make it blank to exclude error rate calc)
				#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
				for($i=0; $i < @stringDSSP; $i++){
					if($stringDSSP[$i] =~ /\W/){ $seq_positionSTR[$i]= $gap_char;}
				}
			}
		}elsif( $comm_col =~ /C/i){
				print __LINE__, " Replacing position with \gap_char \"$gap_char\"\n" if $debug2 eq 1;
				$ss_opt = 'ss'; # whether it was set or not, make it 'ss'
				for($i=0; $i < @stringDSSP_common; $i++){
					if($stringDSSP_common[$i] =~ /\W/){ $seq_positionSTR[$i]= $gap_char;}
				}
		}

		if($debug2 eq 1){
			print __LINE__,
			print " \@seq_positionSTR is  @seq_positionSTR\n";
		}

		#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		#   getting Position differences.
		#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		@position_diffs  = @{&get_posi_diff(\@seq_positionSEQ, \@seq_positionSTR)};

		if($debug2 eq 1){
			print __LINE__,
			print " \@position_diffs is  @position_diffs\n";
		}

		#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		#  You can have two types of output according to which alignment you compare your
		#   error rates. (1) Compare to @stringSEQ   (2) @stringSTR
		#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		@position_corrected1 = @{&put_position_back_to_str_seq(\@stringSEQ, \@position_diffs)};
		$array3{$name}=join(",", @position_corrected1);

	}
	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	# The final Step for error rate, $LIMIT is to confine error rate in one digit (ie, less than 10)
	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	%final_posi_diffs =%{&get_residue_error_rate(\%array3, $LIMIT)};

	undef(@whole_length, $len_of_seq);
	return(\%final_posi_diffs);
}

#________________________________________________________________________
# Title     : convert_num_0_or_1_hash_opposite
# Usage     : with a variable for threshold ->
#
#               %out = %{&convert_num_0_or_1_hash_opposite(\%input_hash, \$threshold)};
#
# Function  : changes all the numbers into 0 or 1 according to threshold given.
#             convert_num_0_or_1_hash converts threshold and bigger nums. to
#             '0' while convert_num_0_or_1_hash_opposite converts to '1'.
# Example   : A hash =>  name1  10012924729874924792742749748374297
#                        name2  10012924729874924792710012924729874
#             A threshold => 4
#             !! if numbers are smaller than 4, they become 1 (or true).
#             Outputhash  =>  name1  11111011011111011111011011110101111
#                        name2  11111011010001011001011010010101100
#
#             ($ref1, $ref2)=&convert_num_to_0_or_1_hash(\%hash, \%hash, \$threshold);
#             above is the example when with more than 2 input hashes.
# Warning   : Threshold value is set to 0 as well as all values smaller than that.
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   :
# Tips      :
# Argument  : two references, one for hash one for scaler for threshold
#
# Todo      :
# Author    : A Biomatic
# Version   : 1.0
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub convert_num_0_or_1_hash_opposite{
  my(@output_hash_refs, %input,$c, $i, $split_char,
	  @string, $name, @names, $threshold,%output_hash);
  for($c=0; $c < @_; $c++){
	 if(ref($_[$c]) eq 'SCALAR'){ $threshold =${$_[$c]};}
	 elsif( $_[$c] =~/^\d+$/){ $threshold = $_[$c];}
  }
  for($i=0; $i <=$#_; $i++){
	 if(ref($_[$i]) eq 'HASH'){
		%input=%{$_[$i]};
		#show_hash(\%input);
		@names=keys (%input);
		$split_char=',';
		if ((@_ == 1)&&(ref($_[$a]) eq 'HASH')){  # if input argument is only one (= if no threshold given),
		  $threshold = 1; } # <---- put 1 to $threshold as a default
		for $name (@names){
		  if($input{$name}=~/\,/){  $split_char = ',';
		  }else{ $split_char = ','; }
		  if ($input{$name} =~ /[\.\-\d]+/){
			 @string=split(/$split_char/, $input{$name});
			 for (@string){
				if(/\d+/){
				  if($_ >= $threshold){ $_ = 1; } # !! becomes 0 (or false)
				  else{  $_=0;               } # !! becomes 1 (or true)
				}
			 }
		  }
		  $output_hash{$name}=join(",", @string);
		}
		push(@output_hash_refs, \%output_hash);
	 }
  }
  if(@output_hash_refs == 1){return($output_hash_refs[0]); }
  elsif(@output_hash_refs > 1){ return(@output_hash_refs) }
}

#________________________________________________________________________
# Title     : parse_arguments
# Usage     : &parse_arguments; or  (file1, file2)=@{&parse_arguments};
# Function  : Parse and assign any types of arguments on prompt in UNIX to
#             the various variables inside of the running program.
#             This is more visual than getopt and easier.
#             just change the option table_example below for your own variable
#             setttings. This program reads itself and parse the arguments
#             according to the setting you made in this subroutine or
#             option table in anywhere in the program.
# Example   : &parse_arguments(1);
#             @files=@{&parse_arguments(1)};
# Warning   : HASH and ARRAY mustn't be like = (1, 2,3) or (1,2 ,3)
# Class     :
# Keywords  :
# Options   : '0'  to specify that there is no argument to sub, use
#              &parse_arguments(0);
#             parse_arguments itself does not have any specific option.
#             '#' at prompt will make a var  $debug set to 1. This is to
#              print out all the print lines to make debugging easier.
# Package   :
# Reference :
# Returns   : Filenames in a reference of array
#             and input files in an array (file1, file2)=@{&parse_arguments};
# Tips      :
# Argument  : uses @ARGV
# Todo      :
# Author    : A Biomatic
# Version   : 1.6
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub parse_arguments{
  my( $c, $d, $f, $arg_num, $option_table_seen, $n, $option_filtered,
		$option_table_example, $input_line, @input_files);
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #   Checks if there were arguments
  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  if( @ARGV < 1 ){ #<-- If Argument is not given at prompt
	  for(@_){
		 if($_ eq '0'){
			 last;
		 }else{
			 print "\n \"$0\" requires at least one Argument, suiciding.\n\n";
			 print chr(7); #<-- This is beeping
			 print "  To get help type \"$0  h\"\n\n\n ";
			 exit;
		 }
	  }
  }
  #""""""""""""""""""""""""""""""""""""""""""""""""""
  #   Some DEFAULT $debug variables for debugging purposes
  #""""""""""""""""""""""""""""""""""""""""""""""""""
  &set_debug_option;
  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #  If there is only one prompt arg. and is [-]*[hH][elp]*, it calls
  #   &default_help and exits
  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  if( ( @ARGV == 1 ) && ($ARGV[0] =~ /^[\-]*[hH\?][elp ]*$/) ){
		&default_help;   exit;
  }
  for($f=0; $f < @ARGV; $f++){
	  if( ($ARGV[$f] =~ /\w+[\-\.\w]+$/)&&(-f $ARGV[$f]) ){
		 push(@input_files, $ARGV[$f] ); next;  }
  }
  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #     Reading the running program script
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""
  &assign_options_to_variables;
  if($HELP == 1){ &default_help }
  return(\@input_files);
}

#________________________________________________________________________
# Title     : print_seq_in_block
# Usage     : &print_seq_in_block (\$block_leng, 'i',\%h1, 'sort', \%h2, \%hash3,,,);
# Function  : gets many refs  for one scalar  or hashes and prints
#               the contents in lines of \$block_leng(the only scalar ref. given) char.
# Options   : 'o' or 'O' => ordered hash print,
#             'n' or'N' => no space between blocks.
#             's' or 'S' => printout sorted by seq names.
#             'i' or 'I' => interlaced print.(this requires identical names in hashes)
#             'v' or 'V' => show sequence start number at each line
#             (all options can be like \$sort
#             while $sort has 's' as value. naked number like 100 will be the
#             block_length. 'i' or 'I' => interlaced print.(this requires
#             identical names in hashes)
# Warning   :
# Class     :
# Keywords  :
# Example   : If there are 3 hashes output will be; (in the order of \%hash3, \%hash2, \%hash1)
#             >> 1st Hash        >> 2nd Hash         >> 3rd Hash
#             Name1  THIS-IS-    Name123  eHHHHHHH   Name123  12222223
#
#             You will get;
#                            Name1    THIS-IS-
#                            Name123  eHHHHHHH
#                            Name123  12222223
#
# Package   :
# Reference :
# Returns   :
# Tips      :
# Argument  : many refs  for hash (one for bottm, one for top, etc,top hash is usually
#               to denote certain caculations or results of the bottom one
# Todo      :
# Author    : A Biomatic
# Version   : 1.1
# Used in   :
# Enclosed  :   -- Following are examples.
#             Example of ( no option, DEFAULT )  # Example of ('i' or 'I' option,
#                                                                INTERLACE )
#             6taa           ----ATPADWRSQSIY    #   6taa       ------ATPADWRSQSIY
#             2aaa           ------LSAASWRTQS    #   6taa       ------CCHHHHCCCCEE
#             1cdg           APDTSVSNKQNFSTDV    #   6taa       ------563640130000
#
#             6taa           ------CCHHHHCCCC    #   2aaa       ------LSAASWRTQSIY
#             2aaa           ------CCHHHHCCCC    #   2aaa       ------CCHHHHCCCCEE
#             1cdg           CCCCCCCCCCCCCCCC    #   2aaa       ------271760131000
#
#             6taa           ------5636401300    #   1cdg       APDTSVSNKQNFSTDVIY
#             2aaa           ------2717601310    #   1cdg       CCCCCCCCCCCCCCCCEE
#             1cdg           6752327236000000    #   1cdg       675232723600000000
#
#             Example of('s' or 'S' option,SORT) # Example of ('o' or 'O' option,
#                                                        ORDERED by input hashes )
#
#             1cdg           APDTSVSNKQNFSTDV    #   6taa       ------ATPADWRSQSIY
#             2aaa           ------LSAASWRTQS    #   2aaa       ------LSAASWRTQSIY
#             6taa           ------ATPADWRSQS    #   1cdg       APDTSVSNKQNFSTDVIY
#
#             1cdg           CCCCCCCCCCCCCCCC    #   6taa       ------CCHHHHCCCCEE
#             2aaa           ------CCHHHHCCCC    #   2aaa       ------CCHHHHCCCCEE
#             6taa           ------CCHHHHCCCC    #   1cdg       CCCCCCCCCCCCCCCCEE
#
#             1cdg           6752327236000000    #   6taa       ------563640130000
#             2aaa           ------2717601310    #   2aaa       ------271760131000
#             6taa           ------5636401300    #   1cdg       675232723600000000
#--------------------------------------------------------------------
sub print_seq_in_block{
  my($c, $d, $e, $f, $k, $i, $s, $t,@in,$gapped,
	  %input0,%input1,%input2,%input3,%input4, @names0, @names1, @names2, @names3,
	  $intl,$z,$diff,$offset);
  my($bl)=60;
  my($sort) =0; my($n_space)=0; my($ordered) =0; my($gap_char) ='-';
  my($na,$larg,$names,$seq, $visual_num); my($n)=13;
  sub num{ $a <=> $b; } my(@in_ar, $bl_passed);

  ###############  ARGV handling ######################
  for($k=0; $k< @_ ;$k++){
	  if( !ref($_[$k]) ){    # when inputs are not ref.
		  if($_[$k]=~ /^[\-]*([\d]{1,4})$/){ $bl=$1 if $1>0; $bl_passed=1 if $1>0; next;}  #<--   option handling
		  if($_[$k]=~ /^[\-sS]$/)   { $sort =    1; next;}
		  if($_[$k]=~ /^[\-nN]$/)   { $n_space = 1; next;}
		  if($_[$k]=~ /^[\-iI]$/)   { $intl =    1; next;}
		  if($_[$k]=~ /^[\-]*[gG]+/){ $gapped   =  1; next;}
		  if($_[$k]=~ /^[\-]*[vV]+/){ $visual_num =  1; next;}
		  elsif($_[$k]=~ /^[\-oO]$/){ $ordered = 1;      }
	  }
	  elsif( ref($_[$k]) eq "SCALAR" ){     #<--   option handling
		  if(${$_[$k]}=~ /^[\-]*([\d]{1,4})$/){$bl=$1 if $1>0; $bl_passed=1 if $1>0; next;}      # the scalar input
		  if(${$_[$k]}=~ /^[\-sS]$/){$sort = 1;next;}                 # should shorter than 5
		  if(${$_[$k]}=~ /^[\-nN]$/){$n_space = 1;next;}              # to be recognized as
		  if(${$_[$k]}=~ /^[\-iI]$/){$intl = 1;next;}                 # options, be it number or
		  if(${$_[$k]}=~ /^[\-]*[g]/i){ $gapped =   1; next;}                 # options, be it number or
		  if(${$_[$k]}=~ /^[\-]*[v]/i){ $visual_num = 1; next;}                 # options, be it number or
		  elsif(${$_[$k]}=~ /^[o]$/i){$ordered = 1;}                  # or chars.
	  }
	  elsif(ref($_[$k]) eq "HASH") {  push(@in,  $_[$k]); } #<-- seqn handling hash
	  elsif(ref($_[$k]) eq "ARRAY"){  push(@in, &convert_arr_and_str_2_hash($_[$k], $k));} #<-- conv array to hash.
	  elsif(ref($_[$k]) eq "SCALAR"){ push(@in,&convert_arr_and_str_2_hash($_[$k], $k));} #<-- conv array to hash.
  }

  #########  HASH input handling ############
  for($k=0; $k< @in; $k++){
		 if(ref($in[$k]) eq "HASH"){
			  %{"input$k"}=%{$in[$k]};
			  print %input0;
			  if($sort == 1){   ## When the keys should be sorted.
				  $keys_long= join("", keys(%{"input$k"}) );   ## makes a string of keys to do the next
				  if ($keys_long =~ /[\d\.]+/){                ## see if there is digit.
					  @{"names$k"}= sort num keys(%{"input$k"}); # numerical sort of keys(seq names)
				  }
				  elsif($keys_long =~ /[\w\.\,]+/){          ## if there is no digits,
					  @{"names$k"}= sort keys(%{"input$k"});  ## do the normal string sort
				  }
			  }elsif($sort != 1){                           ## no sorting at all
					  @{"names$k"}= keys(%{"input$k"});
			  }

			  if($gapped != 1){
				  for($i=0; $i< @{"names$k"}; $i++){
					 if(${"input$k"}{${"names$k"}[$i]} =~ /\,/){               # remove ','
						${"input$k"}{${"names$k"}[$i]}=~ s/\,//g;
					 }
				  }
			  }
		 }
  }

  ########################################################################
  ##     Following is to make ends of sequences neat                   ##
  ########################################################################
  for($z=0; $z < @in; $z++){
	 for($t=0;$t < @{"names$z"}; $t ++ ){
		 $na=${"names$z"}[$t];
		 $s=${"input$z"}{$na};
		 $larg=length($s) if length($s)> $larg;
		 $n=length($na) if length($na) > $n;
		 if($s =~ /\-/){ $gap_char='-'; }elsif( $s =~ /\./){  $gap_char='.';  }
		 if (length($s)<$larg){
			$offset=length($s);
			$diff=$larg-$offset;
			substr($s,$offset,$larg)="$gap_char"x$diff;
		 }
	 }
  }

  ########################################################################
  ##     Following is the core code for making block printing          ##
  ########################################################################
  if($ordered== 1){
		$bl=$larg if (($larg < 60)&&($bl_passed != 1));
		for($c=0; $c < @in; $c++){
			for($k=0; $k < $larg; $k += $bl){
				for($i=0; $i < @{"names$c"}; $i++){
					 $names= ${"names$c"}[$i];
					 $seq= substr( ${"input$c"}{$names},  $k,  $bl);
					 $seq=join(" ", split(/\,/, $seq)) if $gapped == 1;
					 if($visual_num==1){
						 printf ("%-${n}s %4d %-$bl s\n", $names, ($k+1), $seq);
					 }else{ printf ("%-${n}s %-$bl s\n", $names, $seq); }
				}print "\n" unless($n_space == 1);
			}print "\n";
		}
  }elsif($intl==1){   ## When Interlace option is set
		  $bl=$larg  if (($larg < 50)&&($bl_passed != 1));
		  for($k=0; $k < $larg; $k+=$bl){
			 for($i=0; $i < @{"names0"}; $i++){
				for($c=0; $c <= $#in; $c++){
					 $names=${"names$c"}[$i];
					 $seq=substr(${"input$c"}{$names}, $k, $bl);
					 $seq=join(" ", split(/\,/, $seq)) if $gapped == 1;
					 if($visual_num==1){
						 printf ("%-${n}s %4d %-$bl s\n", $names, ($k+1), $seq);
					 }else{ printf ("%-${n}s %-$bl s\n", $names, $seq); }
				}print "\n" unless($n_space ==1);
			 }print "\n";
		  }print "\n" unless($n_space ==1);
  }else{
		############################################################
		##           This is the default                          ##
		############################################################
		for($k=0; $k < $larg; $k+=$bl){
			$bl=$larg if (($larg < 50)&&($bl_passed != 1));
			for($c=0; $c < @in; $c++){  # $n is the name space size
				my(@seq_names) = @{"names$c"};
				for($i=0; $i < @seq_names; $i++){
					 my($names)=$seq_names[$i];
					 my($long_seq)=${"input$c"}{$names};
					 my($seq)=substr($long_seq, $k, $bl);
					 $seq=join(" ", split(/\,/, $seq)) if $gapped == 1;
					 if($visual_num==1){
						 printf ("%-${n}s %4d %-$bl s\n", $names, ($k+1), $seq);
					 }else{ printf ("%-${n}s %-$bl s\n", $names, $seq); }
				}print "\n" unless($n_space ==1);
			}print "\n";
		}
	}
}
#________________________________________________________________________
# Title     : open_dssp_files
# Usage     : (*out, *out2) = @{&open_dssp_files(\$inputfile1, \$inputfile2, \$H, \$S,,,,)};
#             (@out)        = @{&open_dssp_files(\$inputfile1, \$inputfile2, \$H, \$S,,,,)};
# Function  : open dssp files and put sequences in a hash(s)
#              It can take options for specific secondary structure types. For example,
#              if you put an option $H in the args of the sub with the value of 'H'
#              open_dssp_files will only read secondary structure whenever it sees 'H'
#              in xxx.dssp file ignoring any other sec. str. types.
#              If you combine the options of 'H' and 'E', you can get only Helix and long
#              beta strand sections defined as segments. This is handy to get sec. str. segments
#              from any dssp files to compare with pdb files etc.
#             With 'simplify' option, you can convert only all the 'T', 'G' and 'I' sec. to
#              'H' and 'E'.
# Example   :
# Warning   : 6taa.dssp  and 6taa are regarded as the same.
# Class     :
# Keywords  :
# Options   : H, S, E, T, I, G, B, P, C, -help
# $H        =        'H' by   -H or -h or H or h  # to retrieve 4-helix (alpha helical)
# $S        becomes  'S' by   -S or -s or S or s  # to retrieve Extended strand, participates in B-ladder
# $E        becomes  'E' by   -E or -e or E or e  # to retrieve residue in isolated Beta-bridge
# $T        becomes  'T' by   -T or -t or T or t  # to retrieve H-bonded turn
# $I        becomes  'I' by   -I or -i or I or i  # to retrieve 5-helix (Pi helical) segment output
# $G        becomes  'G' by   -G or -g or G or g  # to retrieve 3-helix (3-10 helical)
# $B        becomes  'B' by   -B or -b or B or b  # to retrieve only B segment
# $simplify becomes   1  by   -p or P or -P, p
# $comm_col becomes  'c' by   -c or c or C or -C or common
# $HELP     becomes   1  by   -help   # for showing help
#
# Package   :
# Reference :
# Returns   : (*out, *out2)  or (@out_array_of_refs)
# Tips      :
# Argument  : files names like (6taa, 6taa.dssp) If you put just '6taa' without extension, it
#             searches if there is a '6taa.dssp' in both PWD and $DSSP env. set directory.
#             ---------- Example of dssp ---
#             **** SECONDARY STRUCTURE DEFINITION BY THE PROGRAM DSSP, VERSION JUL
#             REFERENCE W
#             HEADER    RIBOSOME-INACTIVATING PROTEIN           01-JUL-94   1MRG
#             COMPND    ALPHA-MOMORCHARIN COMPLEXED WITH ADENINE
#             SOURCE    BITTER GOURD (CUCURBITACEAE MOMORDICA CHARANTIA) SEEDS
#             AUTHOR    Q
#             246  1  0  0  0 TOTAL NUMBER OF RESIDUES, NUMBER OF CHAINS, NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN)                .
#             112 95.0   ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)                                                                         .
#             171 69.5   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)  , SAME NUMBER PER 100 RESIDUES                              .
#             12   4.9   TOTAL NUMBER OF HYDROGEN BONDS IN     PARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .
#             36  14.6   TOTAL NUMBER OF HYDROGEN BONDS IN ANTIPARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .
#             1    0.4   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-5), SAME NUMBER PER 100 RESIDUES                              .
#             1    0.4   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-4), SAME NUMBER PER 100 RESIDUES                              .
#             74  30.1   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+4), SAME NUMBER PER 100 RESIDUES                              .
#             5    2.0   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+5), SAME NUMBER PER 100 RESIDUES                              .
#             1    2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30     *** HISTOGRAMS OF ***           .
#             0    0  0  0  1  1  0  2  0  0  1  0  0  1  0  0  0  0  0  2  0  0  0  0  0  0  0  0  0  0    RESIDUES PER ALPHA HELIX         .
#             1    0  0  2  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    PARALLEL BRIDGES PER LADDER      .
#             2    0  1  2  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    ANTIPARALLEL BRIDGES PER LADDER  .
#             2    0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    LADDERS PER SHEET                .
#             #   RESIDUE AA STRUCTURE BP1 BP2  ACC   N-H-->O  O-->H-N  N-H-->O  O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
#             1    1   D              0   0  132    0, 0.0   2,-0.3   0, 0.0  49,-0.2   0.000 360.0 360.0 360.0 153.4   44.0   96.9  -23.8
#             2    2   V  E     -a   50   0A  10   47,-1.5  49,-2.8   2, 0.0   2,-0.3  -0.889 360.0-163.3-115.9 151.4   43.1  100.4  -22.5
#             3    3   S  E     -a   51   0A  63   -2,-0.3   2,-0.3  47,-0.2  49,-0.2  -0.961  10.3-172.8-131.0 152.3   44.8  103.7  -23.4
#             4    4   F  E     -a   52   0A   8   47,-2.2  49,-2.3  -2,-0.3   2,-0.4  -0.985   6.9-161.2-143.2 139.5   45.0  107.2  -22.0
#             5    5   R  E     -a   53   0A 144   -2,-0.3   4,-0.2  47,-0.2  49,-0.2  -0.993   9.7-156.0-121.0 125.9   46.6  110.2  -23.6
#             6    6   L  S    S+     0   0    1   47,-2.3   2,-0.5  -2,-0.4   3,-0.4   0.644  73.2  90.9 -73.3 -22.4   47.5  113.2  -21.4
#             7    7   S  S    S+     0   0   81   47,-0.3   3,-0.1   1,-0.2  -2,-0.1  -0.695 106.0   5.2 -75.5 121.0   47.4  115.6  -24.4
#             8    8   G  S    S+     0   0   72   -2,-0.5  -1,-0.2   1,-0.3   5,-0.1   0.269  97.6 147.8  90.2 -10.7   43.9  117.0  -24.7
#             9    9   A        +     0   0   10   -3,-0.4  -1,-0.3  -4,-0.2  -3,-0.1  -0.256  16.8 166.8 -58.8 142.4   42.9  115.2  -21.5
#             (\$inputfile1, \$inputfile2, .... )};
# Todo      :
# Author    : A Biomatic
# Version   : 2.9
#             $debug feature has been added to make it produce error messages with '#' option.
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub open_dssp_files{
	#"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
	my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
	my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
	my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
	if($debug==1){print "\n\t\@hash=\"@hash\"
	\@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
	\@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  my $gap_char = '_';

  if($char_opt !~ /[HEBGIST]/i){  ## This is default sec. str type setting. (full representation)
		$char_opt = 'HEBGIST';
  }
  if ($debug eq 1){
	  print __LINE__, " # open_dssp_files : \$simplify     is  $simplify\n" ;
	  print __LINE__, " # open_dssp_files : \@file given   is  @file \n" ;
	  print __LINE__, " # open_dssp_files : \@string given is  @string\n" ;
  }
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  ### Big main loop for input argument handling   ####
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #  This is to check if the given file is not in pwd but in ENV var $DSSP
  #  Or if the file name was given only by the base name of seq(eg. 1cdg rather
  #  than 1cdg.dssp
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  for ($i=0; $i < @string; $i ++){
		 print __LINE__, " ${i}th string input is  $string[$i] \n" if $debug eq 1;
		 $string[$i] = "$string[$i]\.dssp"; ## adding  .dssp extension
		 print __LINE__, " ${i}th string inputwith \.dssp is now, $string[$i] \n" if $debug eq 1;

		 if(-f $string[$i]){
			 print chr(7) if $debug eq 1;
			 print __LINE__, " Your input filename exist in this File: $string[$i]\n" if $debug eq 1;
			 unshift(@file, "$string[$i]");
		 }
		 elsif(-l $string[$i]){
			 print chr(7) if $debug eq 1;
			 print "\n Your input filename exist as a Link to : $string[$i]\n" if $debug eq 1;
			 unshift(@file, "$string[$i]");
		 }
		 elsif( -d $ENV{'DSSP'} ){
			 $string[$i] =~ s/(\w+)\.dssp$/$1/; ## stripping .dssp extension
			 if( -e "$ENV{'DSSP'}\/$string[$i]\.dssp" ){
				unshift(@file, "$ENV{'DSSP'}\/$string[$i]\.dssp");
				$BASE = $string[$i];
			 }else{
				 print chr(7);
				 print __LINE__, " !! Error your DSSP env setting seems wrong. \n";
				 print __LINE__, " !! Your DSSP env path is also a link. \n" if (-l $ENV{'DSSP'});
				 print __LINE__, " I can't find  $ENV{'DSSP'}\/$string[$i] \n\n";
			 }
		 }
		 elsif( -l $ENV{'DSSP'} ){ #"""""""  IF $DSSP was a link
			 print __LINE__, " !! Your DSSP env path is also a link. \n" if $debug eq 1;
			 if( -e "$ENV{'DSSP'}\/$string[$i]\.dssp" ){
				unshift(@file, "$ENV{'DSSP'}\/$string[$i]\.dssp");
				$BASE = $string[$i];
			 }
			 elsif( -e "$ENV{'DSSP'}\/$string[$i]" ){
				unshift(@file, "$ENV{'DSSP'}\/$string[$i]\.dssp");
			 }
		 }
  }

  @file=@{remove_dup_in_array(\@file)};

  if ($debug eq 1){
	  print __LINE__, " # open_dssp_files : ENV set for dssp is $ENV{'DSSP'} \n" ;
	  print __LINE__, " # open_dssp_files : Final \@file given are \" @file \"\n" ;
  }
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #  END of File and string input checking in searching for the right dssp file.
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #"""""""""""""""""""""""" MAIN """"""""""""""""""""""""""""""""""""""
  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  for($i=0; $i< @file; $i++){  ## <<-- loops over the input files.
		 my($flag, %hash, $name, $s, $matched, $ori_name, $chain);
		 my($real_file) = $file[$i];
		 $file[$i] =~ s/(.*\/)(\w+)\.(\w+)$/$2/; ## stripping .dssp extension
		 $file[$i] =~ s/(\w+)\.(\w+)$/$1/;       ## stripping .dssp extension
		 $ori_name = $name = $file[$i];
		 print "\n",__LINE__, " VAR \$ori_name is  $ori_name , \$file\[\$i\] is $file[$i]\n" if $debug eq 1;
		 unless(-e $real_file){
			print "\n",__LINE__,"  !!! ERROR $real_file does not exists as the final filename\n" if $debug eq 1;
			splice(@file, $i, 1); $i--;
			print "\n",__LINE__,"  Skipping to the next file to open" if $debug eq 1;
			next;
		 }

		 open(FILE_1,"$real_file");
		 print "\n",__LINE__, " ${i}th file $real_file is being opened from \@file \n" if $debug eq 1;
		 print "_"x86,"\n", if $debug eq 1;

		 while(<FILE_1>){
			 if(/^[\s]*\#\s+RESIDUE/){
				 $flag =1;
				 print __LINE__," \"#  RESIDUE\"   string found at line $. in $real_file\n" if $debug eq 1;
				 next
			 }  ##   '#  RESIDUE' is the starting key

			 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
			 #    Matching the column
			 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

			 if(  ($flag==1) && (/^[\s]*-*\d+\s+-*\d*\s+[\w]\s\s([\w ]) /)  ){
				 $matched = $1;
				 print __LINE__," \"$matched\" is matched\n" if $debug2 eq 1;

				 if( $char_opt =~ /$matched/){ ## Here OPTIONS affect the operation.
					 $s .= $matched;    ## $match_option is like 'HE'. If the
					 next;              ## single char $matched is H or E, it will be
				 }else{                ## annexed to $s as an output.
					 $s .= $gap_char;
					 next;
				 }  # <-- this is necessary to get the right length (not to ignore
			 }     #     not matched char by converting them to ' '.

			 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
			 #    When there are chains like A, B, ,,,
			 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
			 elsif( ($flag==1) && (/^[\s]*-*\d+\s+-*(\d+)[\s\w]+(\w)\s+[\w]\s\s([\w ]) /) ){
				 $chain = $2;   ## $flag  is for the starting key
				 # ${"chain_start$name$2"} = $1 unless defined(${"chain_start$name$2"});
				 my($matched_chain) = $3;
				 if( $char_opt =~ /$matched_chain/){
					$s .= $matched_chain;   next; }
				 else{
					$s .= $gap_char; next; }
			 }elsif( (/^\s+\d+\s+\!/)&&($chain =~/\w/) ){
				 $name="$name$chain";
				 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
				 ##   IF simplify  option is set
				 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
				 if($simplify eq 1){
					 $s =~ tr/TGI/EHH/;   ### change the characters.
					 print __LINE__," Simplifying TGI to EHH by \"tr\"\n" if $debug eq 1;
				 }
				 if($debug eq 1){ print __LINE__, " Name of seq:  $name \n"; }
				 $hash{$name}=$s; $s='';
				 $name=$ori_name; next;
			 }
		 }
		 close(FILE_1);  ##<<---- Reading finished.

		 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		 ##  Naming procedure
		 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		 if($chain =~/^\w$/){     # when there are chains, put A,B, etc to seq. names.
			 $name="$name$chain";  ## <<-- This is for the last chain entry.
			 if($debug eq 1){ print __LINE__, " Name of seq:  $name, There were Chains !\n"; }
			 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
			 ##   IF simplify  option is set
			 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
			 if($simplify eq 1){
				 $s =~ tr/TGI/EHH/;   ### change the characters.
				 print __LINE__," Simplifying TGI to EHH by \"tr\"\n" if $debug eq 1;
			 }
			 $hash{$name}=$s;
			 $s='';   ##<<--- This is essential, a former bug!
			 $name=$ori_name;
		 }else{      # <<-- Without chains option.
			 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
			 ##   IF simplify  option is set
			 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
			 if($simplify eq 1){
				 $s =~ tr/TGI/EHH/;   ### change the characters.
				 print __LINE__," Simplifying TGI to EHH by \"tr\"\n" if $debug eq 1;
			 }
			 $hash{$name}=$s;
			 if($debug eq 1){ print __LINE__, " Name of seq:  $name \n"; }
			 $s='';   #<<--- This is essential, a former bug!
		 }

		 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		 ### OUTput format determination according to options #####
		 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
		 if($debug eq 1){ print "\n", __LINE__, " The Hash out of \"$real_file\" is \n ";
			 &show_hash(%hash);
		 }
		 push(@out_hash_ref_list, \%hash) if ref(\%hash) eq 'HASH';
  }
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #"""""""""""" END of Main """""""""""""""""""""""""""
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

  if($comm_col =~ /c/i){ # $comm_col  is a global
	  if($debug eq 1){
		  print "\n", __LINE__;
		  print " # open_dssp_files : you have put 'c' option for common column only\n";
		  $temp = @out_hash_ref_list;
		  print __LINE__, " # open_dssp_files : No. of hashes passed to get_common_column is: $temp\n";
		  print __LINE__, " # open_dssp_files : The hash are(is) : @out_hash_ref_list\n";
	  }
	  $ref_hash_out = &get_common_column(@out_hash_ref_list);
	  return($ref_hash_out);
  }else{
	  if(@out_hash_ref_list == 1){ return($out_hash_ref_list[0]); }
	  elsif(@out_hash_ref_list > 1){ return(@out_hash_ref_list);  }
  }
}

#________________________________________________________________________
# Title     : get_wrong_segment_rate
# Usage     : print_seq_in_block( &get_wrong_segment_rate(\%superposed_hash) );
# Function  : Treats the segment as one single big error.
#             calculates the wrong segment number compared to the correct ones.
# Example   : <input example> hash of 3 keys and values.
#             2aaa_6taa      -------00000---------00000000----0000-------00000-
#             1cdg_6taa      -------442---------------2222-----------------000-
#             1cdg_2aaa      -------222---------------2222-----------------000-
#
#             In the above there are two segments wrong in 3 segment blocks = 2/3
#             <output example> hash of 3 percentage rates.
#
#             2aaa_6taa      0 %
#             1cdg_6taa      66.6666666666667 %
#             1cdg_2aaa      66.6666666666667 %
#
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   : 1.0
# Used in   : get_segment_shift_rate
# Enclosed  :
#--------------------------------------------------------------------
sub get_wrong_segment_rate{
  my($a, $b, $c, $d, $e, $f, $g, $h, $i, $j, $k, $l, $m, $n, $o, $p, $q, $r,
	  $s, $t, $u, $v, $w, $x, $y, $z, %h, $seg_min,
	  %hash, @keys, @array, @hash, $option_string, $string,
	  $name, %out, $gap_chr, @str1, @str2, $seg, $len, $wrong_seg, $correct_seg
  );
  %hash=%{$_[0]};
  $seg_min =$_[1];
  if($seg_min !~/\d+/){ $seg_min = 3; } ### Default segmin is 3
  @keys = sort keys (%hash);
  for $k (@keys){
	 my($string) = $hash{$k}; $string =~s/\,//g;
	 my(@segments) = split(/[\-\.\ ]+/, $string);
	 for $seg (@segments){
		$len=length($seg);
		if( $len >= $seg_min){
			if($seg =~/[1-9]/){
				$wrong_seg ++;  }
			else{ $correct_seg ++; }
		}
	 }
	 $h{$k}= ($wrong_seg/($wrong_seg + $correct_seg)*100).' %';
	 $wrong_seg=$correct_seg='';
  }
  \%h;
}
#________________________________________________________________________
# Title     : get_correct_percent_alignment_rate
# Usage     : &get_correct_percent_alignment_rate(\$file1, \$file2);
# Function  : accepts two files and prints out the sequence identities of the alignment.
# Example   :
# Warning   : Alpha version,  A Biomatic , made for Bissan
# Class     :
# Keywords  :
# Options   : h  # for help
#             v  # for verbose printouts(prints actual sequences)
# Package   :
# Reference :
# Returns   : reference of Scalar for percentage correct alignment(for already
#             aligned sequences)
# Tips      :
# Argument  : two sequence files which have identical sequence names.
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub get_correct_percent_alignment_rate{
	 my($i, $j, $k, $verbose, @string1, @string2, $larger, $seq_pair_id, @seq_pair_ids );
	 my(%inputhash1) = %{&read_any_seq_files($_[0])};
	 my(%inputhash2) = %{&read_any_seq_files($_[1])};
	 my(@names)= sort keys %inputhash1;
	 ######################################
	 ####### Sub argument handling ########
	 ######################################
	 for($k=0; $k< @_ ;$k++){
		if( !ref($_[$k])&& (length(${$_[$k]}) < 5)){  # when inputs are not ref.
		  if($_[$k]=~ /^[\-vV]$/){ $verbose = 1; next;}
		}
		elsif((ref($_[$k]) eq "SCALAR")&&(length(${$_[$k]})<5)){  #  when inputs are  ref.
		  if(${$_[$k]}=~ /^[\-vV]$/){$verbose = 1;next;}          # should shorter than 5
		}
	 }
	 for($i =0; $i < @names; $i++){
		print "\n\n==== Processing structural $names[$i] against artificial $names[$i]\n";
		$inputhash1{$names[$i]} =~ tr/a-z/A-Z/;
		$inputhash2{$names[$i]} =~ tr/a-z/A-Z/;
		@string1=split(//, $inputhash1{$names[$i]});
		@string2=split(//, $inputhash2{$names[$i]});
		print "\n The string1 is  ",@string1,"\n" if $verbose ==1;
		print "\n The string2 is  ",@string2,"\n" if $verbose ==1;
		(@string2 > @string1) ? ($larger=@string2, $smaller=@string1) : ($larger=@string1, $smaller=@string2);
		$true_seq=$inputhash1{$names[$i]};
		$true_seq=~s/\W//g;
		$true_len=length($true_seq);
		print "\n True seq length:   $true_len  , Whole length inc gap: $larger\n";
		for($j = 0; $j < $larger; $j++){
		  $iden_sum++ if ($string1[$j] eq $string2[$j])&& !($string1[$j]=~/\W/); }
		$seq_pair_id =($iden_sum/$true_len) * 100;
		print "\nID between structural and artifical alignment is  $seq_pair_id \%" , "\n";
		push(@seq_pair_ids, $seq_pair_id);
		undef( $iden_sum, $seq_pair_id );
	 }
	 print "\n", "_"x88, "\n";
	 my($whole_average_of_the_id)=${&array_average(\@seq_pair_ids)};
	 print "The whole average is; $whole_average_of_the_id\n";
	 if(@seq_pair_ids == 1){ return( \$seq_pair_ids[0] ); }
	 elsif(@seq_pair_ids > 1){ return( \@seq_pair_ids ); }
}
#________________________________________________________________________
# Title     : read_any_seq_files
# Usage     : %out_seq=%{&read_any_seq_files(\$input_file_name)};
# Function  : Tries to find given input regardless it is full pathname, with or
#             without extension. If not in pwd, it searches the dirs exhaustively.
# Example   : (*out1,  *out2) =&read_any_seq_files(\$input1, \$input2);
#             : (@out_ref_array)=@{&read_any_seq_files(\$input1, \$input2)};
#             : (%one_hash_out) =%{&read_any_seq_files(\$input1)};
# Warning   :
# Class     :
# Keywords  : open_any_seq_files
# Options   : v for $verbose setting showing some information in runtime
#
# Returns   : 1 ref. for a HASH of sequence ONLY if there was one hash input
#             1 array (not REF.) of references for multiple hashes.
# Tips      :
# Argument  : one of more ref. for scalar.
# Todo      :
# Author    : A Biomatic
# Version   : 1.1
#--------------------------------------------------------------------
sub read_any_seq_files{
  my(@out_hash_ref_list, $sub, $o, $ext );
  my(@in)=@_;
  for($o=0; $o< @in; $o++){
	 my($found, %out, @file_ext_accepted, $found_file, $sub);
	 if(ref($_[$o])){
		 @file_ext_accepted=('msf', 'fasta','jp','aln','ali','pir',
								  'slx', 'dna','fas','pdb','rms','brk', 'dssp');
		 if( ! -e ${$in[$o]}  or -B ${$in[$o]} or -z ${$in[$o]}  ){
			 print "\n#SUB read_any_seq_files: ${$in[$o]} no seq file exists(or not passed at all) for $0 \n\n",
			 chr(7);
			 exit;
		 }
		 $found_file=${&find_seq_files($in[$o])};
		 print "# in read_any_seq_files, \$found_file => $found_file\n" if $verbose==1;

		 for $ext(@file_ext_accepted){
			$sub ="open\_$ext\_files";
			print "# Trying subroutine $sub\n" if $verbose==1;
			if($found_file =~/\.$ext$/){
			   %out=%{&{"$sub"}(\$found_file)} if (defined &{"$sub"}); $found =1;
			}
			if($found_file =~/\.$ext$/ and  ! defined &{"$sub"} ){
			   print "\n# $sub is not defined in $0. I want it!!\n\n";
			}
		 }
		 if($found==0){
		    my($sub)="open\_$ext\_files"; #<--- this is the last resort !!
			for $ext(@file_ext_accepted){
			   %out=(%out, %{&{"$sub"}(\$found_file)}) if (defined &{"$sub"});
			}
		 }
	  }elsif( !(ref($_[$o])) ){
	     print "\nread_any_seq_files in $0 files accepts only REFERENCES\n\n";
	     exit;
	  }
	  push(@out_hash_ref_list, \%out);
  }
  if(@out_hash_ref_list == 1){  ### If only single hash output is,
	  return($out_hash_ref_list[0]);
  }elsif( @out_hash_ref_list > 1){
	  return(@out_hash_ref_list);   # <-- contains (\%out_seq0, \%out_seq1, \%out_seq2, .... )
  }
}


#________________________________________________________________________
# Title     : open_jp_files  (bug free!!)
# Usage     : %out_hash=%{&open_jp_files(\$file_name)};
# Function  : reads jp files and stores results in a hash.
# Example   :
# Warning   : All the spaces  '-' !!!
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : a reference of a hash for names and  their sequences.
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   : 1.1
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub open_jp_files{
	#"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
	my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
	my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
	my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
	if($debug==1){print "\n\t\@hash=\"@hash\"
	\@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
	\@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	my(%hash_out, $s1);
	open(FILE_JP, "$file[0]");
	while(<FILE_JP>){    if(/^CLUSTAL/){ next; }
		if((/^([\S]+)[\t]* +$/)||(/^\#/)){ next; }
		if(/^([\w\.\-\=\+]+) +\t*(\S+)[\n]$/){ $n=$1; $s=$2; $hash_out{$n}.= $s; }
	}
	#&show_hash(%hash_out);
	\%hash_out;
}


#________________________________________________________________________
# Title     : open_msf_files
# Usage     : (*out, *out2) = @{&open_msf_files(\$inputfile1, \$inputfile2)};
#             : %hash_seq = %{&open_msf_files(\$inputfile1)};
#             : (@out)        = @{&open_msf_files(\$inputfile1, \$inputfile2)};
#             ---------- Example of MSF ---
#             PileUp
#
#             MSF:   85  Type: P    Check:  5063   ..
#
# Function  : open msf files and put sequences in a hash(s)
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : (*out, *out2)  or (@out_array_of_refs)
# Tips      :
# Argument  : (\$inputfile1, \$inputfile2, .... )};
# Todo      :
# Author    : A Biomatic
# Version   : 1.1
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub open_msf_files{
	#"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
	my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
	my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
	my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
	if($debug==1){print "\n\t\@hash=\"@hash\"
	\@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
	\@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

  my(@names,  %hash, @out_hash_ref_list);
  for($i=0; $i< @file; $i++){
	 open(FILE_1, "$file[$i]");
	 undef(%hash);
	 while(<FILE_1>){      # file1 needs to be xxxx.msf
		 if((/^([\S]+)\t* +$/)||(/^\#/)||(/^\-+/)){ next; }
		 if(/^([\S]+)\t* +([\.\w ]+)[\n]$/){
			 $n=$1;
			 $s=$2;
			 $s=~s/ //g;
			 $hash{$n}.= $s;
		 }
	 }
	 push(@out_hash_ref_list, \%hash);
  }
  if(@out_hash_ref_list  == 1 ){ return(\%hash); }
  elsif(@out_hash_ref_list > 1){ return(@out_hash_ref_list); }
}


#________________________________________________________________________
# Title     : default_help
# Usage     : &default_help2;  usually with 'parse_arguments' sub.
# Function  : Prints usage information and others when invoked. You need to have
#             sections like this explanation box in your perl code. When invoked,
#             default_help routine reads the running perl code (SELF READING) and
#             displays what you have typed in this box.
#             After one entry names like # Function :, the following lines without
#             entry name (like this very line) are attached to the previous entry.
#             In this example, to # Function : entry.
# Example   : &default_help2; &default_help2(\$arg_num_limit);   &default_help2( '3' );
#             1 scalar digit for the minimum number of arg (optional),
#             or its ref. If this defined, it will produce exit the program
#             telling the minimum arguments.
# Warning   : this uses format and references
# Class     :
# Keywords  :
# Options   :
# Package   : File_Util
# Reference :
# Returns   : formated information
# Tips      : This usually goes with  parse_arguments.pl (= easy_opt.pl)
# Argument  :
# Todo      :
# Author    :
# Version   : 3.2
# Used in   : parse_arguments,
# Enclosed  :
#--------------------------------------------------------------------
sub default_help{
  my($i, $perl_dir, $arg_num_limit, $head ,$arg_num_limit );
  my($logname)=getlogin();
  my($pwd)=`pwd`;
  my($date)=`date`;
  chomp($date,$pwd);
  my($not_provided)="--- not provided ---\n";
  my($file_to_read) = $0;

  for($i=0; $i < @_; $i ++){
	  if((ref($_[$i]) eq 'SCALAR')&&(${$_[$i]} =~ /^\d$/)){
		  $arg_num_limit = ${$_[$i]};  }
	  elsif( (!(ref($_[$i]))) && ($_[$i] =~ /^\d$/)){
		  $arg_num_limit = $_[$i];     }
  }
  my %entries = %{&read_head_box(\$file_to_read )};
  if($option_tb_found ==1){
	 @option_tb=@{&read_option_table(\$file_to_read)};
  }
  foreach $help_item (keys %entries){
	  ${$help_item}= $not_provided if( (${$help_item}=~/^[\W]*$/)||( !defined(${$help_item})) );
  }
  #""""""""""""""""""""""""""""""""""""""""
  #########  Writing the format <<<<<<<<<<<
  #""""""""""""""""""""""""""""""""""""""""
  $~ =HEADER_HELP;
  write;   ## <<--  $~ is the selection operator
  $~ =DEFAULT_HELP_FORM;
  for(sort keys %entries){  write  }
  print chr(7);  print "_"x72,"\n\n";

  if(@ARGV < $arg_num_limit){ print "\* $0 fataly needs $arg_num_limit arguments\n\n" }

  if(  $option_tb_found == 1){
	 #########  Printing the OPTION table contents <<<<<<<<<<<<
	 print "  Press \"Return\" key to see what options $logname\'s \n\n    \"$0\" take... \n";
		 $key_press=getc();
	 print @option_tb, "\n"x2 if(@option_tb > 0);
  }
format HEADER_HELP  =
_____________________________________________________________________
		  __  __      ______     __          _____
		 /\ \/\ \    /\  ___\   /\ \        /\  _ `\
		 \ \ \_\ \   \ \ \__/   \ \ \       \ \ \L\ \
		  \ \  _  \   \ \  _\    \ \ \       \ \ ,__/
		   \ \ \ \ \   \ \ \/___  \ \ \_____  \ \ \/
		    \ \_\ \_\   \ \_____\  \ \______\  \ \_\
		     \/_/\/_/    \/_____/   \/______/   \/_/ V 3.1`
_____________________________________________________________________
.
format DEFAULT_HELP_FORM =
 @<<<<<<<<<: @*
 $_        $entries{$_}
.
}
#________________________________________________________________________
# Title     : set_debug_option
# Usage     : &set_debug_option;
# Function  : If you put '#' or  '##' at the prompt of any program which uses
#             this sub you will get verbose printouts for the program if the program
#             has a lot of comments.
# Example   : set_debug_option #    <-- at prompt.
# Warning   :
# Class     : Utility
# Keywords  :
# Options   : #   for 1st level of verbose printouts
#             ##  for even more verbose printouts
# $debug  becomes 1 by '#'  or '_'
# $debug2 becomes 1 by '##'  or '__'
#
# Package   :
# Reference : http://sonja.acad.cai.cam.ac.uk/perl_for_bio.html
# Returns   :  $debug
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   : 1.8
#             generalized debug var is added for more verbose printouts.
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub set_debug_option{
  my($j, $i, $level);
  unless( defined($debug) ){
	 for($j=0; $j < @ARGV; $j ++){
		 if( $ARGV[$j] =~/^(_+)$|^(#+)$/){ # in bash, '#' is a special var, so use '_'
			 print __LINE__," >>>>>>> Debug option is set by $1 <<<<<<<<<\n";
			 $debug=1;
				  print chr(7);
			 print __LINE__," \$debug  is set to ", $debug, "\n";
			 splice(@ARGV,$j,1); $j-- ;
			 $level = length($1)+1;
			 for($i=0; $i < $level; $i++){
				 ${"debug$i"}=1;
				 print __LINE__," \$debug${i} is set to ", ${"debug$i"}, "\n";
			 }
		 }
	 }
  }
}
#________________________________________________________________________
# Title     : remov_com_column
# Usage     : %new_string = %{&remov_com_column(\%hashinput)};
#             @out=@{&remov_com_column(\@array3)};
# Function  : removes common gap column in seq.
# Example   :
# Warning   :
# Class     :
# Keywords  : remove_com_column, remove_common_column,
#             remove_common_gap_column, remov_common_gap_column,
#             remove com column
# Options   :
# Package   :
# Reference :
# Returns   : a ref. of  hash(es) and array(s).
#
#             name1   ABCDE....DDD       name1  ABCDE..DDD
#             name2   ABCDEE..DD..  -->  name2  ABCDEEDD..
#             name3   ...DEE..DDE.       name3  ...DEEDDE.
#
#             (ABC....CD, ABCD...EE) --> (ABC.CD, ABCDEE)
#             from above the two column of dot will be removed
#             To remove absurd gaps in multiple sequence alignment. for nt6-hmm.pl
# Tips      :
# Argument  : accepts reference for hash(es) and array(s).
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub remov_com_column{
  my(@hash_ref_out, $d, $gap_char);
  for($d=0; $d < @_; $d++){
	  if(ref($_[$d]) eq 'HASH'){
	      my($len,@string,@new_string,@string2);
		  my(%input)=%{$_[$d]};
		  my(@common);
		  for (keys %input){
			  @string = split('', $input{$_});
			  if(!(defined(@common))){ @common = (@string);  }
			  else{ for ($k=0; $k < @string; $k++){
				 if (($string[$k] =~ /\W/ )&&($common[$k] =~ /(\W)/)){ $common[$k]= $1;}
				 elsif(($string[$k] =~ /(\W)/)&&(!(defined($common[$k])))){ $common[$k]=$1;}
				 else{ $common[$k]='X';} } } }
		  for (keys %input){ @string2 = split(//, $input{$_});
			  for ($i=0; $i < @string2; $i++){
				 if ($common[$i] ne $string2[$i]){ push(@new_string, $string2[$i]); } }
			  $new_string{$_}= join('', @new_string); @new_string = ();      }
		  push(@hash_ref_out, \%new_string);
	  }
	  elsif(ref($_[$d]) eq 'ARRAY'){
	      my( $k, $y, $x,@string_array, @string);
		  my(@input)=@{$_[$d]};  @common=();
		  for($y=0; $y< @input; $y++){
			  @string = split('', $input[$y]);
			  if(!(defined(@common))){  @common = @string;  }
			  else{
				 for ($k=0; $k < @string; $k++){
					 if (($string[$k]  =~ /(\W)/)&&($common[$k]  =~ /(\W)/)){ $common[$k]=$1;}
					 elsif(($string[$k] =~ /(\W)/)&&(!(defined($common[$k])))){ $common[$k]=$1;}
					 else{ $common[$k]='X';}
				 }
			  }
		  }
		  for($x=0; $x < @input; $x++){
		      my($new_string, @string2);
			  @string2 = split(//, $input[$x]);
			  for ($i=0; $i < @string2; $i++){
				  if ($common[$i] ne $string2[$i]){ $new_string.= "$string2[$i]"; }
			  }
			  push(@string_array, $new_string);
		  }
		  push(@hash_ref_out, \@string_array);
	  }
  }
  if(@hash_ref_out ==1) { return( $hash_ref_out[0] ); }
  elsif(@hash_ref_out>1){ return( @hash_ref_out ) }
}
#________________________________________________________________________
# Title     : tidy_secondary_structure_segments
# Usage     : print_seq_in_block(&tidy_secondary_structure_segments(\%hash, 'e4', 'h4'), 's');
#
# Function  : receives any secondary structure assignment hashes and
#             tidys up them. That is removes very shoft secondary structure
#             regions like( --HH--, -E-, -EE- ) according to the given minimum
#             lengths(threshold) of segments by you.
# Example   : print_seq_in_block(&tidy_secondary_structure_segments(\%hash, 'e4', 'h4'), 's');
#             <makes following into the next block>
#
#             1cdg_2aaa      -------EEE-----------EE--EEEE------EE---------EEE-
#             1cdg_6taa      -------EEE-----------EE--EEEE------EE---------EEE-
#             2aaa_6taa      -------EEEEE------EE-EEEEEEEE----EEEE-------EEEEE-
#
#             <example output>
#
#             1cdg_6taa      -------------------------EEEE---------------------
#             1cdg_2aaa      -------------------------EEEE---------------------
#             2aaa_6taa      -------EEEEE---------EEEEEEEE----EEEE-------EEEEE-
#
# Warning   :
# Class     :
# Keywords  :
# Options   : something like 'H3' or 'E3' for minimum segment length set to 3 positions.
# Package   : Bio::Seq
# Reference :
# Returns   : array of references of hashes.
# Tips      :
# Argument  : hashes and [options]. No options result in default of 'H3', 'E3'
# Todo      :
# Author    : A Biomatic
# Version   : 1.0.0
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub tidy_secondary_structure_segments{
  my($i, $k,$a, $j, $helix_min, $beta_strand_min, %hash, @keys, @hash,
	  $option_string, @hash_out, $string1, $name, %out, $gap_chr, @str1, @str2,
	  @stringout, @string_segH, @string_segE, $countH, $countE
	  );

  #### Default helix and beta strand segment length setting #####
  $helix_min=3;
  $beta_strand_min=3;

  ########################################################################
  #####   general argument handling  for options of segment length  ######
  ########################################################################
  for($k=0; $k< @_ ;$k++){
	  if( ( !ref($_[$k]) )&&($_[$k]=~ /^[Hh](\d+)$/) ){
		  $helix_min  = $1;    }
	  elsif( ( !ref($_[$k]) )&&($_[$k]=~ /^[Ee](\d+)$/) ){
		  $beta_strand_min  = $1;    }
	  elsif((ref($_[$k]) eq "SCALAR")&&(${$_[$k]}=~ /^[Hh](\d+)$/) ){
		  $helix_min  = $1;    }
	  elsif((ref($_[$k]) eq "SCALAR")&&(${$_[$k]}=~ /^[EeBb](\d+)$/) ){
		  $beta_strand_min  = $1;    }
	  elsif(ref($_[$k]) eq "HASH") { push(@hash,  $_[$k]); }    }

  for($i=0; $i < @hash; $i++){
	  my(%hash) = %{$hash[$i]};
	  @keys = sort keys( %hash );
	  for($j=0; $j < @keys; $j++){
		  my(@string_segH, @string_segE, @stringout);
		  $string1=$hash{$keys[$j]};
		  $gap_char = $1 if ($string1=~ /(\W)/);

		  ##### actual cleaning ####
		  my(@string) = split(//, $string1);
		  for($a = 0; $a < @string; $a++){
			 if($string[$a] !~/[HE]/){ ### if the splited element doesn't match 'H' or 'E'

				 ##### If any of the HH or EE counter is over the given minimum($helix_min,,)
				 if((@string_segH >= $helix_min)||( @string_segE >=$beta_strand_min)){
					 push(@stringout, @string_segH, @string_segE, '-');
					 @string_segH=@string_segE=();     }   ## just resetting.
				 else{  ### if the accumulated 'HH' or 'EE' is smaller than the minimum
					 for(0.. (@string_segH + @string_segE) ){
						push(@stringout, '-'); ### replace the short 'EE' etc with '-'
					 }
					 @string_segH=@string_segE=();  ## just resetting.
				 }
			 }
			 elsif($string[$a] =~ /^([Hh])$/){
				 push(@string_segH, $1); }
			 elsif($string[$a] =~ /^([Ee])$/){
				 push(@string_segE, $1); }
		  }
		  $hash{$keys[$j]}=join("", @stringout);
	  }
	  push(@hash_out, \%hash);
  }
  if(@hash_out == 1){ return($hash_out[0]);
  }elsif(  @hash_out > 1 ){ return(@hash_out); }
}
#________________________________________________________________________
# Title     : get_posi_diff
# Usage     : @position_diffs =&get_posi_diff(\@seq_position1,\@seq_position2);
# Function  :
# Example   : @compacted_posi_dif =(1 ,2, 1, 1, '.' ,2,  1,  1, '.');
#             @compacted_posi_dif2=(4 ,2, 1, 1, ,2,  1, '.' ,3,  1);
#             output ==> ( 3 0 0 0 . 1 . 2 .)   (it ignores positions which have non digits.
#             output ==> (-3 0 0 0 . 1 .-2 .) when abs is not used.
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : one ref. for an @array of differences of input arrays. array context.
# Tips      :
# Argument  : Takes two ref. for arrays which have positions of residues.
# Todo      :
# Author    : A Biomatic
# Version   : 1.4
# Used in   : evalign.pl, get_position_shift_rate
# Enclosed  :
#--------------------------------------------------------------------
sub get_posi_diff{
	my(@positions1)=@{$_[0]};
	my(@positions2)=@{$_[1]};
	my(@num_diffs_between_str_and_ali, $diff, $z, $gap_char);
	if($debug eq 1){
	  print __LINE__, " # get_posi_diff : \n";
	}
	$gap_char = '.';
	for ($z=0; $z < @positions2; $z++){
	  if (($positions1[$z] =~ /\d+/) && ($positions2[$z] =~ /\d+/)){
		  $diff=($positions1[$z] - $positions2[$z]);
		  push(@num_diffs_between_str_and_ali, $diff );
	  }else{
		  push(@num_diffs_between_str_and_ali, $gap_char);
	  }
	}
	\@num_diffs_between_str_and_ali;
}
#________________________________________________________________________
# Title     : get_posi_sans_gaps
# Usage     : @seq_position1 = &get_posi_sans_gaps($string1);
# Function  :
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : the positions of residues after removing gaps(but keeps pos).
#               used for analysis of shifted positions of seq. comparison.
# Tips      :
# Argument  : one scalar variable input of sequence string.
# Todo      :
# Author    : A Biomatic
# Version   : 1
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub get_posi_sans_gaps{
  my($string) = ${$_[0]};
  my($char, @positions, $i);
  for($i=0; $i < length($string); $i++){
	 $char=substr($string,$i,1);
	 if(($char eq '-')||($char eq '.')){  next; }else{ push(@positions, $i); } }
  \@positions;
}
#________________________________________________________________________
# Title     : get_common_column
# Usage     : %out =%{&get_common_column(\%hash1, \%hash2, '-')};
# Function  : (name1         --EHH--HHEE-- )
#             (name2         --HHH--EEEE-- ) ==> result is;
#
#             (name1_name2   -- HH--  EE-- )
#             to get the identical chars in hash strings of sequences.
#
# Example   : %out =%{&get_common_column(\%hash1, \%hash2, '-')};
#             output> with 'E' option >>> "name1     --HHH--1232-"
#   Following input will give;
#   %hash1 = ('s1', '--EHH-CHHEE----EHH--HHEE----EHH--HHEE----EHH-CHHEE--');
#   %hash2 = ('s2', '--EEH-CHHEE----EEH-CHHEE----EEH-CHHEE----EEH-CHHEE--');
#   %hash3 = ('s3', '-KEEH-CHHEE-XX-EEH-CHHEE----EEH-CHHEE----EEH-CHHEE--');
#   %hash4 = ('s4', '-TESH-CHEEE-XX-EEH-CHHEE----EEH-CHHEE----EEH-CHHEE--');
#
#     s1_s2_s3_s4    --E-H-CH-EE----E-H--HHEE----E-H--HHEE----E-H-CHHEE--
#
# Warning   : This gets more than 2 hashes. Not more than that!
#
# Class     : get_common_column, get_common_column_in_seq,
#             get common column in sequence, superpose_secondary_structure,
#             get_common_secondary_structure,
#             for secondary structure only representation.
# Keywords  : Overlap, superpose hash, overlay identical chars, superpose_seq_hash
#             get_common_column, get_com_column, get_common_sequence,
#             get_common_seq_region, multiply_seq_hash,
# Options   :
# Package   : Array_Util
# Reference :
# Returns   : one hash ref. of the combined key name (i.e., name1_name2). Combined by '_'
# Tips      :
# Argument  : 2 or more ref for hash of identical keys and value length.
#             One optional arg for replacing space char to the given one.
# Todo      :
# Author    : A Biomatic
# Version   : 1.5
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub get_common_column{
  my($i, $k,$j, $name1, $name2, @in, %out, @out_chars, $gap_chr, @str1, @str2);
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #  Sub argument handling     $gap_chr here can be 'HE' etc.
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""
  for($k=0; $k< @_ ;$k++){
	  if( ( !ref($_[$k]) )&&($_[$k]=~ /^(.)$/) ){
		  $gap_chr  .= $1;    }
	  elsif((ref($_[$k]) eq "SCALAR")&&(${$_[$k]}=~ /^(.)$/) ){
		  $gap_chr  .= $1;    }
	  elsif(ref($_[$k]) eq "HASH") { push(@in,  $_[$k]); }    }

  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #"  Checking if the input hashes were right
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""
  if( (@in < 2) && ( ref($in[0]) eq 'HASH') ){
	  print "\n", __LINE__, " # get_common_column usually needs 2 hashes. Error \n";
	  print "\n", __LINE__, " # get_common_column : Anyway, I will just return the single input hash:  @in. \n";
	  %out=%{$in[0]}; # <--- This is essential to return the single input hash!!
	  goto FINAL;
  }

  %out = %{$in[0]};  ## Initializing %out
  print "\n",__LINE__, " # get_common_column hashes given are: @in \n" if $debug eq 1;

  for( $k=1; $k < @in; $k++){
		my(@out_chars);   ## <-- Necessary to prevent concatenation.
		my(%hash1)=%out;
		my(%hash2)=%{$in[$k]};
		my(@names1)= sort keys %hash1;
		my(@names2)= sort keys %hash2;
		$name1 = $names1[0];
		$name2 = $names2[0];
		@str1=split(/||\,/, $hash1{$names1[0]});
		@str2=split(/||\,/, $hash2{$names2[0]});
		for($i=0; $i < @str1; $i++){
			if($str1[$i] eq $str2[$i] ){
				push(@out_chars, $str1[$i]); }
			elsif( defined($gap_chr) ){ push(@out_chars, $gap_chr); }
			else{ push(@out_chars, ' '); }
		}
		if( $name2 < $name1){      ## To make an ordered name output eg.  seq1_seq2, than  seq2_seq1
			%out='';
			$out{"$name2\_$name1"}= join("", @out_chars); }
		else{
			%out='';
			$out{"$name1\_$name2"}= join("", @out_chars); }
  }
  FINAL:
  if ($debug eq 1){
	  print "\n",__LINE__, " # get_common_column Final res. \%out :\n",
	  &show_hash(%out);
  }
  \%out;
}
#________________________________________________________________________
# Title     : array_average
# Usage     : $output = &array_average(\@any_array);
# Function  : (the same as average_array)
# Example   :
# Warning   : If divided by 0, it will automatically replace it with 1
# Class     :
# Keywords  : get_array_average, av_array, average_array, get_average_array
#             average_of_array, average_array
# Options   :
# Package   :
# Reference :
# Returns   : single scaler digit.
# Tips      :
# Argument  : takes one array reference.
# Todo      :
# Author    : A Biomatic
# Version   : 1.2
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub array_average{
  my(@input)= @{$_[0]};
  my $int_option = ${$_[1]} || $_[1];
  my($item,$average,$num,$sum);
  my $num_of_elem = @input;

  for $item(@input){
	 if( $item =~ /^$/ ){  ## If it matches nothing. '$item == 0' does not work !!!
		$num_of_elem --; ## This is to make sure that the denominator does not
	 }                  ## count blank element. (to get correct element number)
	 else{ $sum += $item;  }
  }
  if($num_of_elem ==0){ $num_of_elem =1; }  ## To prevent 'Division by 0' error
  if($int_option =~ /[\-]*i[nt]*/){
	  $average= int( $sum/$num_of_elem );
  }else{   $average = $sum/$num_of_elem }

  return(\$average);
}
#________________________________________________________________________
# Title     : find_seq_files
# Usage     : $found_file = ${&find_seq_files(\$input_file_name)};
# Function  : (similar to find.pl) used in 'read_any_seq_file.pl'
#             seeks given test file in pwd, specified dir, default path etc.
#             If not found yet, it looks at all the subdirectories of path and pwd.
#             PATH environment dirs, then returns full path file name.
# Example   : $found_file=${&find_seq_files(\$input_file_name)};
# Warning   :
# Class     :
# Keywords  : find_anyj_seq_files, find any seq files, find seq files
# Options   :
# Package   :
# Reference :
# Returns   : return( \$final );
# Tips      :
# Argument  : (\$input_file_name) while $input_file_name can be  'xxx.xxx', or '/xxx/xxx/xxx/xxy.yyy'
#             or just directory name like 'aat' for  /nfs/ind4/ccpe1/people/A Biomatic /jpo/align/aat
#             then, it tries to find a file with stored seq file extensions like msf, jp, pir etc
#             to make aat.msf, aat.jp, aat.pir ... and searches for these files.
# Todo      :
# Author    : A Biomatic
# Version   : 1.0
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub find_seq_files{
  my($final, $no_ext_file, $result); my($in_file)=${$_[0]}; my($pwd)=`pwd`; chomp($pwd);
  my( $base, @ENV_dir, $ext, @probable_dir_list, $directory);
  my(@extension_db)=('sst','msf','fasta','jp','fas','aln','brk','pdb', 'rms', 'ent','slx','fa');
  @probable_dir_list=('JPO','ALIGN','PATH','HOME','PIRDIR','PWD','PDBSST','PDBENT','BLASTDB','PIRDIR','SWDIR','PDB');
	if(($in_file=~/\//)&&(-e $in_file)){ $final=$in_file; }
	elsif((-e $in_file)&&(-s $in_file)&&($in_file !~/\//)){ $in_file="$pwd\/$in_file"; $final=$in_file;}
	######## if it was like  '/nfs/ind4/ccpe1/people/A Biomatic /perl.msf'
	elsif($in_file =~ /\/([\w\-\.]+)$/){ $in_file = $1;
		  if(-e $in_file){ $final = "$pwd\/$in_file"; }
		  #### if it has xxxxxx.xxxx  file form. #######
		  elsif($in_file =~ /(([\w\-]+)\.([\w\-]+))$/){ $file=$1; $base=$2; $ext=$3;
				for (@extension_db){ if($_ eq $ext){ shift(@extension_db);}}
				unshift(@extension_db, $ext);
				for (@probable_dir_list){ if($ENV{$_}=~ /\/$/){chop($ENV{$_});}
					push( @ENV_dir, split(/:/, $ENV{$_}));}
					for $dir (@ENV_dir){ $in_file="$dir\/$file";
						if ((-e $in_file) && (-s $in_file)){  $final=$in_file; last;}
						else{
							 for $ext (@extension_db){ $in_file="$dir\/$base\.$ext";
								  if ((-e $in_file) && (-s $in_file)){
									  if ($file =~  /$in_file/){ $final = $in_file; last;}}}}}
					unless(defined ($final)){
						for $dir (@ENV_dir){ $in_file= ${&search_files_in_subdir(\$dir, \$file)};
							 if(-e $in_file){ $final=$in_file; last; }}}}

			### if it has  xxxxxx   file form, ie. not extension #######
			elsif($in_file =~/\/([\w_\-]+)$/){  $base = $1;
			  for (@extension_db){
				 if($_ eq $ext){ shift(@extension_db);  }
				 unshift(@extension_db, $ext);
				 for (@probable_dir_list){
					if ($ENV{$_} =~ /\/$/){  chop($ENV{$_}); }
					push( @ENV_dir, split(/:/, $ENV{$_}) );
					for $dir (@ENV_dir){ $no_ext_file="$dir\/$base";
						 if((-e $no_ext_file) && (-s $no_ext_file)){ $final=$no_ext_file; last;}
						 else{
							for $ext (@extension_db){ $in_file ="$dir\/$base\.$ext";
								if ((-e $in_file) && (-s $in_file)){ $final = $in_file; last;}}}}}}}}

	 #### when the input was like this  'perl.msf'  in any directory.
	 elsif($in_file =~ /^(([\w\-]+)\.([\w\-]+))$/){ $file=$1; $base=$2; $ext=$3;
		  for (@extension_db){ if($_ eq $ext){ shift(@extension_db);}}
		  unshift(@extension_db, $ext);
		  for (@probable_dir_list){ if($ENV{$_}=~ /\/$/){chop($ENV{$_});}
			  push( @ENV_dir, split(/:/, $ENV{$_}));}
			  for $dir (@ENV_dir){ $in_file="$dir\/$file";
				  if ((-e $in_file) && (-s $in_file)){ $final=$in_file; last;}
				  else{
						for $ext (@extension_db){ $in_file="$dir\/$base\.$ext";
							 if ((-e $in_file) && (-s $in_file)){
								 if ($file =~  /$in_file/){ $final = $in_file; last;}}}}}
			  unless(defined ($final)){
				  for $dir (@ENV_dir){ $in_file= ${&search_files_in_subdir(\$dir, \$file)};
						if(-e $in_file){ $final=$in_file; last; }}}}
	 #### when the input was like this  'hemocyan'  in any directory.
	 elsif($in_file =~ /^([\w\-]+)$/){ $file=$1;
		  for (@probable_dir_list){ if($ENV{$_}=~ /\/$/){chop($ENV{$_});}
			  push( @ENV_dir, split(/:/, $ENV{$_}));}
			  for $dir (@ENV_dir){ $in_file="$dir\/$file";
				  if ((-e $in_file) && (-T $in_file)){  $final=$in_file; last;}
				  else{
						for $ext (@extension_db){ $in_file="$dir\/$file\.$ext";
							 if ((-e $in_file) && (-s $in_file)){
								 if ($file =~  /$in_file/){ $final = $in_file; last;}}}}}
			  unless(defined ($final)){
				  for $dir (@ENV_dir){ $in_file= ${&search_files_in_subdir(\$dir, \$file)};
						if(-e $in_file){ $final=$in_file; last; }}}}
	END_POINT:
	return( \$final );
}
#________________________________________________________________________
# Title     : put_position_back_to_str_seq
# Usage     : @result =@{&put_position_back_to_str_seq(\@string_from_struct, \@compacted_posi_dif)};
# Function  :
# Example   : @string_from_struct=('X', 'T', 'A' ,'B' , '.' ,'F',  'G', '.' , 'O' ,'P', '.');
#             @compacted_posi_dif=(1 ,2, 1, 1, ,2, 1, 1, 1);
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : a ref. for an array
# Tips      :
# Argument  : takes two refs for arrays (one for char the other for digits
# Todo      :
# Author    : A Biomatic
# Version   : 1.0
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub put_position_back_to_str_seq{
  my(@string_from_struct)=@{$_[0]};
  my(@compacted_posi_dif)=@{$_[1]};
  my($j)=0; my($char)=0; my($i);
  for ($i=0; $i < @string_from_struct; $i++){
	 $char = $string_from_struct[$i];
	 if ($char =~ /\w/){
		 $string_from_struct[$i] = $compacted_posi_dif[$i-$j];
	 }else{ $j++; }
  }
  \@string_from_struct;
}
#________________________________________________________________________
# Title     : hash_common_by_keys
# Usage     : %hash1_value = %{&hash_common_by_keys(\%hash1, \%hash2,...)};
# Function  :
# Example   :
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : the VALUES OF THE FIRST HASH which occur in later hashes
#             are returned
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub hash_common_by_keys{  my(%common)=();
  for($i=0; $i< @_; $i++){  my(%common2)=();
	  if( !(defined(%common) )){ %common=%{$_[$i]}; next;}
	  elsif(defined(%common)){ %h1=%{$_[$i]};
		 for(keys %common){ $common2{$_}=$common{$_} if (defined $h1{$_});}
	  %common=%common2;}
	  undef(%common2);  }
  \%common;
}
#________________________________________________________________________
# Title     : convert_arr_and_str_2_hash
# Usage     : ($hash1, $hash2)=&convert_arr_and_str_2_hash(\$input, \$input2, '1', '2'.. );
#             * This is the combination of convert_string_to_hash & convert_array_to_hash
# Function  : makes hash(es) out of array(s)
#             if ordering digit(s) is put, it orders the keys according to it.
#             if ordering digit is not increased by one, the difference is used
#             as the increasing factor. No option results in
#             array_1, array_2, array_3...
#
# Example   : &print_seq_in_block(&convert_arr_and_str_2_hash(\@input,\@input2,\@input3 ));
#             &convert_arr_and_str_2_hash(\$input1,\$input2, '2' );
#             results in; (ordering starts from the given '2')
#                          array_2       input1 arraystring
#                          array_3       input2 arraystring
#
#             one more exam
#                          string_6       This is st                  and 3 strings)
#                          string_10      This is st
#                          array_2        111233434242
#                          array_6        111233434242
#                          array_10       111243424224
# Warning   :
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : one or more ref. of hashes.
# Tips      :
# Argument  : one or more ref. of arrays
# Todo      :
# Author    : A Biomatic
# Version   :
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub convert_arr_and_str_2_hash{ my(@out_ref_hash_str,@s,$i,@w,@c,$diff,$c);
  undef(%temp); undef(%temp2);
  for($i=0; $i<@_; $i++){
	 if ( (ref($_[$i]) eq 'SCALAR')&&(${$_[$i]}=~/^([\d]{1,5})$/) ){
		 push(@c, $1);      }
	 elsif( (ref($_[$i]) eq 'SCALAR')){
		 push(@s, ${$_[$i]});     }
	 elsif( (ref($_[$i]) eq 'ARRAY')){
		 push(@w, join("", @{$_[$i]}));    }
	 elsif( (!ref($_[$i])) && ($_[$i]=~/^([\d]{1,5})$/)){
		 push(@c, $1);      }
	 else{ print "\n There is an odd arg. check convert_arr_and_str_2_hash in $0\n\n";
		print chr(007);   exit;  }  }
  for($i=0; $i<@s; $i++){ my($string)=$s[$i]; my(%temp); ###### array handling
	  if(defined($c[$i])){ $temp{"string_$c[$i]"}=$string;
		  push( @out_ref_hash_str, \%temp); $c=$c[$i]; $diff =$c[$i]-$c[$i-1];      }
	  elsif( !(defined($c[$i])) ){
		  if($diff ==0){	$c=$c+$diff+1; }else{ ($c=$c+$diff) };
		  $temp{"string_$c"}=$string; push( @out_ref_hash_str, \%temp);   } }
  for($i=0; $i<@w; $i++){ my($string)=$w[$i]; my(%temp2);###### string handling
	  if(defined($c[$i])){ $temp2{"array_$c[$i]"}=$string; push( @out_ref_hash_str, \%temp2);
		  $c=$c[$i]; $diff =$c[$i]-$c[$i-1];      }
	  elsif( !(defined($c[$i])) ){
		  if($diff ==0){	$c=$c+$diff+1; }else{ ($c=$c+$diff) };
		  $temp2{"array_$c"}=$string;  push( @out_ref_hash_str, \%temp2); } }
  if( @out_ref_hash_str == 1 ){ return($out_ref_hash_str[0] ); }
  elsif(@out_ref_hash_str > 1){ return(@out_ref_hash_str);}
}
#________________________________________________________________________
# Title     : get_residue_error_rate
# Usage     : %position_diffs =%{&get_residue_error_rate(\@seq_position1, \@seq_position2)};
# Function  : This is the final step in error rate getting.
#             gets a ref. of a hash and calculates the absolute position diffs.
# Example   :
# Warning   : split and join char is ',';
# Class     :
# Keywords  :
# Options   : 'L' for limitting the error rate to 9 to make one digit output
#  $LIMIT becomes 'L' by L, l, -l, -L
# Package   :
# Reference :
# Returns   : one ref. for an array of differences of input arrays. array context.
#             ---Example input (a hash with sequences); The values are differences after
#                                comparion with structural and sequential alignments.
#             %diffs =('seq1', '117742433441...000',   <-- input (can be speparated by '' or ','.
#                      'seq2', '12222...99999.8888',
#                      'seq3', '66222...44444.8822',
#                      'seq4', '12262...00666.772.');
#             example output;
#             seq3_seq4       '0,1,0,0,0,.,.,.,,.,0,,0,0,,0,0,,.,0,,0,0,.'
#             seq1_seq2       '0,1,0,1,1,.,.,.,,.,2,,2,2,,2,2,,.,.,,2,2,1'
#             seq1_seq3       '0,1,0,1,1,.,.,.,,.,1,,1,1,,0,.,,.,.,,1,1,1'
#             seq1_seq4       '0,1,0,,1,1,.,.,.,,.,1,,1,1,0,.,.,,.,1,,2,2'
#             seq2_seq3       '0,1,0,,0,0,,.,.,,.,0,,1,0,,0,0,,.,0,,0,0,0'
#             seq2_seq4       '0,0,0,,1,0,,.,.,,.,0,,1,0,,0,0,,.,0,,0,0,.'
# Tips      :
# Argument  : Takes a ref. for hash which have positions of residues of sequences.
# Todo      :
# Author    : A Biomatic
# Version   : 1.1
# Used in   : get_position_shift_rate, previously get_each_posi_diff_hash
# Enclosed  :
#--------------------------------------------------------------------
sub get_residue_error_rate{
	my ($LIMIT);
	my(%diffs)= %{$_[0]}; my(@names)= keys (%diffs);
	$LIMIT=${$_[1]} if ref($_[1]) eq 'SCALAR';
	$LIMIT= $_[1] unless ref($_[1]);
	my(%seqs_comp_in_pair, @temp, @temp2,$split_char, $i);
	for ($i=0; $i < @names; $i++){
		if($diffs{$names[$i]}=~/\,/){ $split_char =',';}else{ $split_char = ''; }
		(@{"string$i"}) = split(/$split_char/, $diffs{$names[$i]});   }
	for ($i=0; $i < @names; $i++){
		for ($j=$i+1; $j < @names; $j ++){
			for ($k=0; $k < @string0; $k++){
				if ((${"string$i"}[$k] =~ /[-\d+]/) && (${"string$j"}[$k] =~ /[-\d+]/)){
					my($diff) = abs(${"string$i"}[$k] - ${"string$j"}[$k]);
					if( ($LIMIT =~/L/i)&&($diff > 9) ){ push(@temp2, 9);
					}else{ push(@temp2, $diff); }
				}else{ push(@temp2, '.'); } }

			#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
			#  Following if {} is for sorting output names to make  2aaa_6taa than 6taa_2aaa
			#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
			if($names[$i] <= $names[$j]){
				$seqs_comp_in_pair{"$names[$i]\_$names[$j]"}=join(",", @temp2); }
			else{ $seqs_comp_in_pair{"$names[$j]\_$names[$i]"}=join(",", @temp2); }

			@temp2=();
		}
	 }
	\%seqs_comp_in_pair;  # permutated
}
#________________________________________________________________________
# Title     : handle_arguments
# Usage     : Just put the whole box delimited by the two '###..' lines below
#             to inside of your subroutines. It will call 'handle_arguments'
#             subroutine and parse all the given input arguments.
#             To use, claim the arguments, just use the variable in the box.
#             For example, if you had passed 2 file names for files existing
#             in your PWD(or if the string looks like this: xxxx.ext),
#             you can claim them by $file[0], $file[1] in
#             your subroutine.
# Function  : Sorts input arguments going into subroutines and returns default
#             arrays of references for various types (file, dir, hash, array,,,,)
#             If you give (\@out, @file), it will put @out into @array as a ref
#             and also the contents of @out will be dereferenced and put to
#             raw_string regardless what is in it).
#
# Example   : 'handle_arguments(\@array, $string, \%hash, 8, 'any_string')
# Warning   :
# Class     : Perl::Utility::Arg_handling
# Keywords  : handling arguments, parsing arguments,
# Options   :
# Package   :
# Reference :
# Returns   : Following GLOBAL variables
#
#             $num_opt,    @num_opt     @file          @dir
#             $char_opt,   @char_opt    %vars          @array,
#             @hash        @string,     @raw_string    @range,
#
#             $num_opt has 10,20
#             @num_opt has (10, 20)
#             @file has  xxxx.ext
#             @dir has  dir  or /my/dir
#             $char_opt has 'A,B'
#             @char_opt has (A, B)
#             @array has  (\@ar1, \@ar2)
#             @hash has (\%hash1, \%hash2)
#             @string  ('sdfasf', 'dfsf')
#             @raw_string (file.ext, dir_name, 'strings',,)
#             @range has values like  10-20
#             %vars deals with x=2, y=3 stuff.
#
# Tips      : takes 0.02 u time with INDY
# Argument  : any type, any amount
# Todo      :
# Author    : A Biomatic
# Version   : 4.6
#             set_debug_option  is added.
# Used in   : everywhere
# Enclosed  :
#--------------------------------------------------------------------
sub handle_arguments{
	my($c, $d, $e, $f, $i, $j, $k, $l, $s, $t, $x, $y, $z, $char_opt, $dir, @hash,
		$file, $in_dir, $num_opt, @char_opt, @dir, @file, @string, @file_dir, @k,
		@num_opt, @raw_string,@string, @array, %vars, @range, @temp, $temp,
		@char_options);

  &set_debug_option;
  if(@_<1){ print chr(7),"\n This is handle_arguments. No args Passed, Error?\n"}
  elsif( (@_ ==1)&& (ref($_[0]) eq 'ARRAY') ){ # when there is only 1 argument
	  push(@array, $_[0]);
	  push(@k, $_[0]);
  }elsif( (@_==1)&&( !ref($_[0]) ) ){
	  if(-f $_[0]){ push(@file, $_[0]);   push(@string, $_[0]) }
	  elsif(-d $_[0]){ push(@dir, $_[0]); push(@string, $_[0]) }
	  elsif($_[0]=~/^\d+$/){ push(@num_opt, $_[0]); $num_opt.=$_[0] }
	  elsif($_[0]=~/^\w+$/){ push(@string, $_[0]); }
  }elsif(@_ >=1){ @k = @_ }

  #####______Start of  general argument handling______######
  for($k=0; $k < @k ;$k++){
	  if( !ref($k[$k]) ){
		  if($k[$k]=~ /^[\-]?([a-zA-Z]\d*) {0,5}$/){  push(@char_opt, $1); $char_opt .= "$1\,";
		  }elsif($k[$k]=~ /^\-([a-zA-Z]+)$/){          ## When multiple option is given,
			  @char_options = split(/\,|/, $1);  push(@char_opt, @char_options);
			  $char_opt .= join("\,", @char_options); ## '-' should be used. eg. '-HEGI'
		  }elsif($k[$k]=~ /^(\w+)\=(\S* *)$/){  $vars{$1}=$2;  $vars .= "$1\,";
		  }elsif($k[$k]=~ /^(\-?\d+)$/){ push(@num_opt, $1);  $num_opt .= "$1\,";
		  }elsif($k[$k]=~ /^\d+\.?\d*\-\d+\.?\d*$/){  push(@range,  $k[$k] );
		  }elsif(-f $k[$k]){                          push(@file,   $k[$k] );
		  }elsif(-d $k[$k]){                          push(@dir,    $k[$k] );
		  }elsif($k[$k]=~ /\/[\w\d\.\-]+[\/].+[\/]/){ push(@dir,    $k[$k] );
		  }elsif($k[$k]=~ /^\/[\w\d\.\-]+[\/]*$/){    push(@dir,    $k[$k] );
		  }elsif($k[$k]=~ /^[\/\w\d\-\.]+\.\w+$/){    push(@file,   $k[$k] );
		  }elsif($k[$k]=~/^\w+[\w\d\.\-]+$/){         push(@string, $k[$k] ); # string does not have space
		  }else{                                      push(@raw_string, $k[$k] );  }

	  }elsif( ref($k[$k]) ){
		  if( ref($k[$k]) eq "SCALAR"){
			 if(${$k[$k]} =~ /^[\-]?([a-zA-Z]\d*) {0,5}$/){ push(@char_opt, $1); $char_opt  .= "$1\,";
				}elsif(${$k[$k]}=~ /^\-([a-zA-Z]+)$/){ push(@char_opt, @char_options);
					$char_opt  .= join("\,", @char_options);  ## as an option string.
				}elsif(${$k[$k]}=~ /^(\w+)\=(\S* *)$/){  $vars{$1}=$2;  $vars .= "$1\,";
				}elsif(${$k[$k]}=~ /^(\-?\d+)$/){ $num_opt .= "$1\,";  push(@num_opt, $1);
			    }elsif(${$k[$k]}=~ /^\d+\.?\d*\-\d+\.?\d*$/){    push(@range,  $k[$k] );
				}elsif(-f ${$k[$k]}){                            push(@file,   ${$k[$k]} );
				}elsif(-d ${$k[$k]}){                            push(@dir,    ${$k[$k]} );
				}elsif(${$k[$k]}=~ /\/[\/\w\d\.\-]+[\/].+[\/]/){ push(@dir,    ${$k[$k]} );
				}elsif(${$k[$k]}=~/^\/[\/\w\d\.\-]+[\/]*$/){     push(@dir,    ${$k[$k]} );
				}elsif(${$k[$k]}=~ /^[\/\w\d\-\.]+\.\w+$/){      push(@file,   ${$k[$k]} );
				}elsif(${$k[$k]}=~/^\w+[\w\d\.\-]+$/){           push(@string, ${$k[$k]} );
				}else{                                           push(@raw_string, ${$k[$k]}); }
		  }elsif(ref($k[$k]) eq "ARRAY"){ my @temp_arr = @{$k[$k]}; push(@array, $k[$k]);
			for ($i=0; $i<@temp_arr; $i++){
			   if(-f $temp_arr[$i]){                            push(@file, $temp_arr[$i]);
			   }elsif($temp_arr[$i]=~/^\d+\.?\d*\-\d+\.?\d*$/){ push(@range,$temp_arr[$i] );
			   }elsif(-d $temp_arr[$i]){                        push(@dir , $temp_arr[$i]);
			   }elsif($temp_arr[$i]=~/\/[\/\w\d\.\-]+[\/].+[\/]/){ push(@dir, $temp_arr[$i] );
			   }elsif($temp_arr[$i]=~/^\/[\/\w\d\.\-]+[\/]*$/){ push(@dir, $temp_arr[$i] );
			   }elsif($temp_arr[$i]=~/^[\/\w\d\-\.]+\.\w+$/){   push(@file,$temp_arr[$i] );
			   }elsif($temp_arr[$i]=~/^\w+[\w\d\.\-]+$/){       push(@string,$temp_arr[$i]);
			   }else{                                           push(@raw_string, $temp_arr[$i]); }
			 }
		  }elsif(ref($k[$k]) eq "HASH"){                             push(@hash,   $k[$k] ); }
	  }
  }
  @raw_string=(@raw_string, @string);
  @file = @{&remove_dup_in_arrayH(\@file)};
  #-----------------------------------------------------
	 sub remove_dup_in_arrayH{  my($i, @out_ref, %duplicate, @orig, @out_ref);
		for($i=0; $i<@_; $i++){  undef(%duplicate);
	 if(ref($_[$i]) eq 'ARRAY'){    @orig = @{$_[$i]};    }
	 @nondup = grep { ! $duplicate{$_}++ } @orig; push(@out_ref, \@nondup);  }
		if(@out_ref ==1){ return($out_ref[0]);}
		elsif(@out_ref >1){  return(@out_ref);}
	 }
  #-----------------------------------------------------
  return(\@hash, \@array, \@string, \@dir, \@file, \@num_opt,
			\@char_opt, \$num_opt, \$char_opt, \@raw_string, \%vars, \@range );
}
#________________________________________________________________________
# Title     : show_hash
# Usage     : &show_hash(\@input_array);
# Function  : for debugging purpose. Shows any array elem line by line.
#             the line is 60 elements long (uses recursion)
# Example   : Output:      item1
#             Output:      item2
#             Output:      item3
# Warning   : There is a global variable:  $show_hash_option
#             It tries to detect any given sting which is joined by ','
# Class     :
# Keywords  :
# Options   : -s or -S or s or S for spaced output. Eg)
#             seq1       1 1 1 1 1 1 1 1 1 1 1 1
#
#             instead of
#             seq1       111111111111
#
#             -h or -H or h or H for horizontal line of '---------...'
#
# Package   : Array_Util
# Reference :
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Biomatic
# Version   : 1.7
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub show_hash{
  my($k, $i, $t, @in2, $in, $LEN, %TEM ); ## You should not put $show_hash_option
  my(@in)=@_;                     ## and $horizontal_line  in my !!!
  my($KL)=2; # default keys string length;
  my($VL)=80; # default values string length;
  my($GAP)=2;  # default space between keys and values
  my($horizontal_line, $show_hash_optionXX, $Hash_counter, @line);

  ## This is to get the option of 'space' to make spaced output.
  for($t=0; $t < @in; $t++){
	 if($in[$t] =~/^[-]+[sS][pace]*$/){
		 $show_hash_optionXX = 1;
		 splice(@in, $t, 1);
	 }elsif(${in[$t]} =~/^[-]+[sS][pace]*$/){
		 $show_hash_optionXX = 1;
		 splice(@in, $t, 1);
	 }elsif($in[$t] =~/^[-]+[hH][rR]*$/){
		 $horizontal_line = 1;
		 splice(@in, $t, 1);
	 }elsif(${in[$t]} =~/^[-]+[hH][rR]*$/){
		 $horizontal_line = 1;
		 splice(@in, $t, 1);
	 }
  }

  ######## Main loop #################
  if($horizontal_line ==1){  ## This puts the delimiter '--------------(  )'
	  $Hash_counter ++;
	  print "\n","-"x78,"(${Hash_counter}th hash)", "\n";
  }

  for($k=0; $k < @in; $k++){
	 if(ref($in[$k]) eq 'ARRAY'){  ### When the hashes were given in array ref.
		 &show_hash(@{$in[$k]}, $show_hash_optionXX, $horizontal_line);
		 print "\n";
	 }
	 elsif(ref($in[$k]) eq 'HASH'){  ### recursion
		 &show_hash(%{$in[$k]});
		 print "\n";
	 }
	 elsif(ref($in[$k+1]) eq 'HASH'){  ### recursion
		 &show_hash(%{$in[$k+1]}); print "\n";
	 }
	 elsif(ref($in[$k]) eq 'SCALAR'){  print ${$_[$k]}, "\n";  }
	 elsif( !ref($in[$k]) ){
		 if( !ref($in[$k+1]) && defined($in[$k+1])  ){
			 if($show_hash_optionXX == 1){  #### space option checking.
				#if($in[$k+1] =~ /\,.+\,/){  #### if the string is joined with ','
				#	 @line = split(/\,/, $_[$k+1]);
				# }else{
				#	 @line = split(//, $_[$k+1]);
				# }
				%TEM = @in;
				$LEN = ${&max_elem_string_array_show_hash(keys %TEM)};
				 if($LEN > $KL){ $KL = $LEN + $GAP +2};
				 printf ("%-${KL}s ", $in[$k]);  $k++;
				 printf ("%-${VL}s\n","@line");
			 }else{                        ### If not option is set, just write
				%TEM = @in;
				$LEN = ${&max_elem_string_array_show_hash( keys %TEM)};
				 if($LEN > $KL){ $KL = $LEN + $GAP +2};
				 printf ("%-${KL}s ", $in[$k]);  $k++; # print $in[$k], "\t";  $k++;
				 printf ("%-${VL}s\n",$in[$k]);        # print $in[$k], "\n";
			 }
		 }
		  #________________________________________________________
		  # Title    : max_elem_string_array_show_hash
		  # Keywords : largest string length of array
		  # Function : gets the largest string length of element of any array of numbers.
		  # Usage    : ($out1, $out2)=@{&max_elem_array(\@array1, \@array2)};
		  #            ($out1)       =${&max_elem_array(\@array1)          };
		  # Argument : numerical arrays
		  # returns  : one or more ref. for scalar numbers.
		  # Version  : 1.1
		  #-------------------------------------------------------
		  sub max_elem_string_array_show_hash{
			 my(@input, $i, $max_elem);
			 @input = @{$_[0]} || @_;
			 for($i=0; $i< @input ; $i++){
					$max_elem = length($input[0]);
					if (length($input[$i]) > $max_elem){
						$max_elem = length($input[$i]);
					}
			 }
			 \$max_elem;
		  }
		  #####################################insert_gaps_in_seq_hash
	 }
  }
}
#________________________________________________________________________
# Title     : remove_dup_in_array
# Usage     : @out2 = @{&remove_dup_in_array(\@input1, \@input2,,,,)};
#             @out1 = @{&remove_dup_in_array(\@input1 )};
# Function  : removes duplicate entries in an array.
# Example   : (1,1,1,1,3,3,3,3,4,4,4,3,3,4,4);  --> (1,3,4);
# Warning   :
# Class     :
# Keywords  : merge array elements, remove_repeting_elements,
#             remove_same_array_elements
# Options   :
# Package   :
# Reference :
# Returns   : one or more references.
# Tips      :
# Argument  : one or more refs for arrays or one array.
# Todo      :
# Author    :
# Version   : 1.3
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub remove_dup_in_array{
  my($i, $sort_opt, @out_ref, @nondup,%duplicate, @orig, @out_ref);
  my @in=@_;
  for($i=0; $i<@in; $i++){
	 if($in[$i] eq 's'){
		$sort_opt=1;
		splice(@in, $i, 1);
		$i--;
	 }elsif( (ref($in[$i]) eq 'SCALAR')&&(${$in[$i]} eq 's') ){
		$sort_opt=1;
		splice(@in, $i, 1);
		$i--;
	 }
  }
  for($i=0; $i<@in; $i++){
	  undef(%duplicate);
	  if(ref($in[$i]) eq 'ARRAY'){    @orig = @{$in[$i]};    }
	  else{ @orig=@in }
	  @nondup = grep { ! $duplicate{$_}++ } @orig;    ## NOTE -> $_
	  if($sort_opt==1){ @nondup= sort @nondup }
	  push(@out_ref, \@nondup);  }
  if(@out_ref ==1){ return($out_ref[0]);}
  elsif(@out_ref >1){  return(@out_ref);}
}
#________________________________________________________________________
# Title     : assign_options_to_variables
# Usage     : &assign_options_to_variables(\$input_line);
# Function  : Assigns the values set in head box to the variables used in
#             the programs according to the values given at prompt.
#             This produces global values.
#             When numbers are given at prompt, they go to @num_opt
#              global variable. %vars global option will be made
#
# Example   : When you want to set 'a' char to a variable called '$dummy' in
#             the program, you put a head box commented line
#             '#  $dummy    becomes  a  by  -a '
#             Then, the parse_arguments and this sub routine will read the head
#             box and assigns 'a' to $dummy IF you put an argument of '-a' in
#             the prompt.
# Warning   : This is a global vars generator!!!
# Class     :
# Keywords  :
# Options   : '#' at prompt will make a var  $debug set to 1. This is to
#              print out all the print lines to make debugging easier.
# Package   : Bio::Utils
# Reference :
# Returns   : Some globaly used variables according to prompt options.
#             @num_opt,
#
# Tips      : Used with 'parse_arguments'
# Argument  : None.
# Todo      :
# Author    : A Scientist
# Version   : 2.4
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub assign_options_to_variables{
  my($i, $j, $op, $z, $n, $symb, $value, $var, %val, @val, $option_table_example, @input_options);

  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #      Defining small variables for option table reading
  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  my($g)='gets';                my($if)='if';
  my($is)='is';                 my(@input_files);
  my($o)='or';   my(@arguments) = sort @ARGV;

  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #  Assigning global arguments(@num_opt, %vars) variables
  #_______________________________________________________________
  for($i=0; $i< @arguments; $i++){
	 if(($arguments[$i]=~/^(\-?\d+[\.\d+]?)$/)&&   ### it mustn't be a file
		( !(-f $arguments[$i]) ) ){                ### getting NUM opt
		push(@num_opt, $1);
	 }elsif( $arguments[$i]=~/^(\S+)=(\S+)$/){
		$vars{$1}=$2;
	 }
  }

  #""""""""""""""""""""""""""""""""""""""""""""""""""
  #   Some DEFAULT $debug variables for debugging purposes
  #""""""""""""""""""""""""""""""""""""""""""""""""""
  &set_debug_option;

  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #   The main processing of self
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  open(SELF, "$0");    ## opens the program you ran to get the options table.
  while(<SELF>){

	  if( $first_border_and_title > 6 ){  ## This is to make it read only the first headbox.
		  last;                            #  $first_border_and_title is an incremental counter.
	  }elsif( (/^ *#[_\*\-]{15,}$/) || (/^ *# *[Tt][itle]*[ :]*/) ){
		  $first_border_and_title++;
		  print __LINE__, "# assign_options_to_variables : Title line found\n" if $debug eq 1;
	  }elsif(/^ {0,5}# {1,50}[\$\%\@].+$/){
		  $op = $&;  ## $op is for the whole input option line which has $xxxx, @xxx, %xxxx format
		  $op =~ s/^( *\# *)(\W\w+.+)$/$2/;  ## This is removing '#  ' in the line.
		  $op =~ s/^(\W\w+.+)(\s+\#.*)$/$1/;  ## This is removing any comments in the line.
			 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
			 ## matching the following line input format.
			 ## $av_sc_segment     becomes    a  by  a  # To smooth the SC rates. Gets the averages of
			 ## $ARG_REG is for arguments regular expression variable.
			 ##  This reg. exp. matches = 'a or A or E or e' part
			 ##  which represents alternative prompt arguments possibilities. \=$b$g$is$e$set
			 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
			 $ARG_REG ='(\S*) *[or=\,]* *(\S*) *[=or\,]* *(\S*) *[=or\,]* *(\S*) *[=or\,]* *(\S*)';
			 if($op=~/^([\$\@\%])([\w\-]+) {0,20}[=|$g|$is] *[\$\@\%]*([\- \w\.\d]+) *[bB]y +$ARG_REG/){
							 ## $sym     $var        becomes          a [$a...]       by       a -a -A
				  my $sym = $1;  #### The symbols like ($, @, %), '$' in the above.
				  my $var = $2;  #### Actual variable name 'var' from $var, 'av_sc_segment' in the above.
				  my $val = $3;  #### The becoming value  first 'a' in the above.
				  my @arg = ($4, $5, $6, $7, $8);  ## The alternative prompt arguments, second 'a' in the above..
			      print "\n $sym $var $val \n" if $debug==1;
			      print "\n \@arg are @arg \n" if $debug==1;

				  #""""""""""""""""""""""""""""""""""""""""""""""""""""
				  #  Going through the PROMPT args.
				  #""""""""""""""""""""""""""""""""""""""""""""""""""""
				  for($z=0; $z < @arguments; $z++){     ## $arguments[$z]  is from @ARGV
					  $arguments[$z] =~ s/\-//;
					  for ($i=0; $i < @arg; $i ++ ){
						 if( ("$arg[$i]" eq "$arguments[$z]" )&& ($arg[$i] !~ /\=/)
							 && ($sym eq '$') ){
							 ${"$var"}="$val";
							 if($debug == 1){
								 print __LINE__," \$${var} is set to \"$1\"\n";
							 }

						 }#'''''''''''''''' $arg = by s=  syntax ~~~~~~~~~~~~~~~~~~~~~~~~~~~
						 elsif( ( $arg[$i] =~ /^(\w+) *\=/ ) &&
							( $arguments[$z] =~ /^${1}= *([\w\.*\-*]+)$/) &&
							( $sym eq '$') ){
							  ${"$var"}="$1";
							  if($debug eq 1){ print __LINE__,"\$${var} is set to \"$1\"\n";  }
						 }
					  }
				  }
			  }
		}
	}
}
