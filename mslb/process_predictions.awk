BEGIN{
  actual = -1
  max_tp = -1
  max_phase = -1
}

function assert(condition, string)
{
    if (! condition) {
        printf("%s:%d: assertion failed: %s\n",
            FILENAME, FNR, string) > "/dev/stderr"
        _assert_exit = 1
        exit 1
    }
}

/Actual object loads/ {
  phase = $5
  if(phase > max_phase){
    max_phase = phase
  }
  actual = 1
}

/Object load predictions/ {
  phase = $5
  if(phase > max_phase){
    max_phase = phase
  }
  actual = 0
}

/tp .* load .*/ {
  tp = $2
  load = $4

  if(tp > max_tp){
    max_tp = tp
  }

  assert((actual >= 0), "Bad phase mutex!")

  #print "actual " actual " tp " tp " phase " phase " load " load

  tp_loads[tp,phase,actual] = load
}

/Done actual object loads/ {
  donePhase = $6
  assert((phase == donePhase), "Bad phase!");
  actual = -1
}

/Done Object load predictions/ {
  donePhase = $6
  assert((phase == donePhase), "Bad phase!");
  actual = -1
}

END{
  for(i=0; i <= max_tp; i++){
    for(j=0; j <= max_phase; j++){ 
      pred_load = tp_loads[i,j,0]
      actual_load = tp_loads[i,j,1]

      print " tp " i " phase " j " predicted " pred_load " actual " actual_load
      delta = (pred_load-actual_load)
      if(delta < 0){
        delta = -delta
      }

      if((actual_load > 0) && (delta/actual_load > err_tol)){
       # print "tp " i " phase " j " predicted " pred_load " actual " actual_load
      }
    }
  }


  if (_assert_exit) exit 1
}
