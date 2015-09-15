# convert.awk : converts out.latslons file to simple format of one 
# ordered pair of depth/soundspeed per line, with profiles separated
# by a line with a ">" character.

BEGIN {
  profile_num=0;
  new_profile=0;
  oldlat=oldlon=-999;
  level=1;
}

substr($1,1,1)!="%" {
  lat=$1; lon=$2; year=$3; month=$4; day=$5; time=$6;
  ssp=$10; depth=$7;

  if( lat!=oldlat || lon!=oldlon || time!=oldtime || day!=oldday || \
    month!=oldmonth || year!=oldyear) {  # then new profile
    profile_num++;
    new_profile=1;
  }

  print olddepth, oldssp;

  if( level>=3 && new_profile && olddepth<=500) {
    # new_profile - eg, just finished last profile so time to mark it.
    # Variables in this section must refer back to previous record,
    # because we don't know we've hit the end of a profile till we've
    # passed it...
    print ">";
    level=1;  # reset level counter and correct first level
  }

  oldlat=lat; oldlon=lon;
  oldtime=time; oldday=day; oldmonth=month; oldyear=year;
  oldssp=ssp;
  olddepth=depth; level++; new_profile=0;
}

