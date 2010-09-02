BEGIN{  
}

/\[edge\] /{
  dist=$5
  vol=$4
  if(vol >= pare){
    print dist " " vol 
  }
}

END{
}
