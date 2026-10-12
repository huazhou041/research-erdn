#!/bin/bash
# Self-healing monitor scoped to the rnaseq_202609 folder ONLY (identifies jobs by
# WorkDir so the leftover source-folder natcommad jobs are never touched).
# Re-kicks node-fault "held/launch failure" jobs via idempotent resubmit, logs
# progress, writes DONE marker when the folder's queue fully drains.
FOLDER=/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
cd "$FOLDER"
LOG=self_heal_202609.log
DONE=DONE_ALL_202609.marker
rm -f "$DONE"
: > "$LOG"
resubmit() {
  bash run_submit_20260709_glm_fsmethods.sh   >>"$LOG" 2>&1
  bash submit_rfxgb.sh                         >>"$LOG" 2>&1
  bash submit_glm_natcommad.sh                 >>"$LOG" 2>&1
  bash submit_glm_smoking.sh                   >>"$LOG" 2>&1
  bash submit_rfxgb_smoking.sh                 >>"$LOG" 2>&1
  bash submit_glm_natcommad_smoking.sh         >>"$LOG" 2>&1
}
while true; do
  ts=$(date '+%F %T')
  # our jobs = those whose WorkDir is this folder
  mapfile -t rows < <(squeue -u zhouh05 -h -o "%i|%t|%r|%Z" 2>/dev/null | awk -F'|' -v d="$FOLDER" '$4==d {print}')
  total=${#rows[@]}
  run=0; pend=0; held=0; heldids=""
  for r in "${rows[@]}"; do
    st=$(echo "$r" | cut -d'|' -f2); rs=$(echo "$r" | cut -d'|' -f3); id=$(echo "$r" | cut -d'|' -f1)
    case "$st" in R) run=$((run+1));; PD) pend=$((pend+1));; esac
    if echo "$rs" | grep -qiE 'launch failure|held|requeue'; then held=$((held+1)); heldids="$heldids $id"; fi
  done
  dres=$(ls res/*_group.done 2>/dev/null | wc -l)
  dsmk=$(ls res_smoking/*_group.done 2>/dev/null | wc -l)
  drf=$(ls res_rfxgb/*_group.done 2>/dev/null | wc -l)
  drfs=$(ls res_rfxgb_smoking/*_group.done 2>/dev/null | wc -l)
  echo "$ts total=$total R=$run PD=$pend held=$held | done: res=$dres res_smoking=$dsmk res_rfxgb=$drf res_rfxgb_smoking=$drfs" >>"$LOG"
  # self-heal: if progress stalled by held jobs, cancel them and refill
  if [ "$held" -gt 0 ] && [ "$run" -eq 0 ]; then
    echo "$ts HEAL: cancelling held [$heldids] and resubmitting" >>"$LOG"
    for id in $heldids; do scancel "$id" 2>>"$LOG"; done
    sleep 20
    resubmit
  fi
  if [ "$total" -eq 0 ]; then
    echo "$ts ALL DONE (folder queue empty) | done: res=$dres res_smoking=$dsmk res_rfxgb=$drf res_rfxgb_smoking=$drfs" >>"$LOG"
    touch "$DONE"
    break
  fi
  sleep 300
done
