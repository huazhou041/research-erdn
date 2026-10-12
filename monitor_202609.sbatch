#!/bin/bash
#SBATCH --partition=cpu_long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --job-name=heal609
#SBATCH --output=self_heal_202609.slurm.%j.log
# Self-healing + SELF-RENEWING monitor scoped to rnaseq_202609 (by WorkDir).
# DONE requires ALL expected group.done counts (empty queue alone is NOT done).
# Renews itself before its 7-day wall so a monitor timeout can never leave the run
# unwatched. A fresh (non-renewed) launch truncates the log + clears markers; a
# renewed launch (RENEWED=1) appends and preserves markers.
FOLDER=/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
cd "$FOLDER"
LOG=self_heal_202609.log
DONE=DONE_ALL_202609.marker
STALL=STALLED_202609.marker
if [ -z "${RENEWED:-}" ]; then rm -f "$DONE" "$STALL"; : > "$LOG"; else echo "$(date '+%F %T') [renewed monitor $SLURM_JOB_ID resumed]" >>"$LOG"; fi
EXP_RES=480; EXP_RFXGB=240; EXP_SMK=384; EXP_RFXGB_SMK=192
MAX_REFILL=6; refill_attempts=0
RENEW_AFTER=594000   # ~6.875 days of script runtime -> renew before 7d wall
resubmit() {
  bash run_submit_20260709_glm_fsmethods.sh >>"$LOG" 2>&1
  bash submit_rfxgb.sh                       >>"$LOG" 2>&1
  bash submit_glm_natcommad.sh               >>"$LOG" 2>&1
  bash submit_glm_smoking.sh                 >>"$LOG" 2>&1
  bash submit_rfxgb_smoking.sh               >>"$LOG" 2>&1
  bash submit_glm_natcommad_smoking.sh       >>"$LOG" 2>&1
}
selfid=$SLURM_JOB_ID
while true; do
  ts=$(date '+%F %T')
  # self-renew before wall-time
  if [ "$SECONDS" -ge "$RENEW_AFTER" ]; then
    echo "$ts RENEW: approaching wall, submitting successor" >>"$LOG"
    sbatch --export=ALL,RENEWED=1 monitor_202609.sbatch >>"$LOG" 2>&1
    break
  fi
  squeue -u zhouh05 -h -o "%i|%t|%r|%Z" 2>/dev/null | awk -F'|' -v d="$FOLDER" -v me="$selfid" '$4==d && $1!=me {print}' > /tmp/heal609.$$
  total=$(wc -l < /tmp/heal609.$$)
  run=$(awk -F'|'  '$2=="R"{c++} END{print c+0}'  /tmp/heal609.$$)
  pend=$(awk -F'|' '$2=="PD"{c++} END{print c+0}' /tmp/heal609.$$)
  held=$(awk -F'|' 'tolower($3) ~ /launch failure|held|requeue/{c++} END{print c+0}' /tmp/heal609.$$)
  heldids=$(awk -F'|' 'tolower($3) ~ /launch failure|held|requeue/{print $1}' /tmp/heal609.$$)
  rm -f /tmp/heal609.$$
  dres=$(ls res/*_group.done 2>/dev/null | wc -l)
  dsmk=$(ls res_smoking/*_group.done 2>/dev/null | wc -l)
  drf=$(ls res_rfxgb/*_group.done 2>/dev/null | wc -l)
  drfs=$(ls res_rfxgb_smoking/*_group.done 2>/dev/null | wc -l)
  echo "$ts total=$total R=$run PD=$pend held=$held | done res=$dres/$EXP_RES res_smk=$dsmk/$EXP_SMK rfxgb=$drf/$EXP_RFXGB rfxgb_smk=$drfs/$EXP_RFXGB_SMK" >>"$LOG"
  if [ "$dres" -ge "$EXP_RES" ] && [ "$drf" -ge "$EXP_RFXGB" ] && [ "$dsmk" -ge "$EXP_SMK" ] && [ "$drfs" -ge "$EXP_RFXGB_SMK" ]; then
    echo "$ts ALL DONE (all expected group.done counts met)" >>"$LOG"; touch "$DONE"; break
  fi
  if [ "$held" -gt 0 ] && [ "$run" -eq 0 ]; then
    echo "$ts HEAL: cancel held [$heldids] + resubmit" >>"$LOG"
    for id in $heldids; do scancel "$id" 2>>"$LOG"; done
    sleep 20; resubmit
  elif [ "$total" -eq 0 ]; then
    refill_attempts=$((refill_attempts+1))
    if [ "$refill_attempts" -le "$MAX_REFILL" ]; then
      echo "$ts REFILL ($refill_attempts/$MAX_REFILL): queue empty but counts short -> resubmitting missing" >>"$LOG"
      resubmit; sleep 60
    else
      echo "$ts STALLED after $MAX_REFILL refills -- surfacing shortfall" >>"$LOG"; touch "$STALL"; break
    fi
  else
    refill_attempts=0
  fi
  sleep 300
done
