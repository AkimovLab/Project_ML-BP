
1. Edit run.py

   * trajectory file name
   * istep, fstep, njobs
   * cp2k template name


2. Edit submit_template.slm

   * number of hours, memory
   * results directory (absolute path)
   * minband, maxband, homo index - should be consistent with the cp2k template
   * the number of states in the "python" script

3. Edit cp2k input template

   * make sure you create enough CUBE files
   * make sure you have enough added MOs
   * make sure the number of excited states is consistent with the one in the submit_template



