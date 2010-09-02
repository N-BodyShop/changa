BEGIN{
  first = 1
}

/\[edge\] /{
  data = $4
  if(first){
    first = 0
    min = data
    max = data
  }
  else{
    if(max < data){
      max = data 
    }
    if(min > data){
      min = data 
    }
  }

}

END{
  print min " " max
}
