BEGIN{
  _gran = (max-min)/nb
  gran = int(_gran)
  if(_gran > gran){
    gran++
  }

  print "min: " min
  print "max: " max
  print "nb: " nb
  print "gran: " gran
  print "pare: " pare

  for(i=0; i<nb; i++){
    counts[i] = 0;
  }
}

/\[edge\] /{
  vol = $4
  if(vol >= pare){
    bin = int((vol-min)/gran)
      counts[bin]++
  }
}

END{
  for(i=0; i<nb; i++){
    imin = min+i*gran
    if(imin > max){
      break
    }
    imax = imin+gran
    if(imax > max){
      imaxstr = max"] "
    }
    else{
      imaxstr = imax") "
    }
    print "[" imin "," imaxstr counts[i]
  }
}
