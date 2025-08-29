new_data_signal()

if not adata:
  w_text_output(content="   ")
  exit()

wf_bigwigs_signal()

if wf_bigwigs_signal.sample() == True:
  tracks = []
  for f in files:
      tracks.append({
          "name": f.path.split("/")[-1],
          "type": "wig",
          "format": "bigwig",
          "url": f.path,
          "autoscale": False,
          "visibilityWindow": 100000,
      })
  
  opts: IGVOptions = {
      "genome": wf_genome.value,
      "tracks": tracks,
  }
  
  igviewer = w_igv(options=opts)

else:
  w_text_output(
    content="   ",
  )
  submit_widget_state()