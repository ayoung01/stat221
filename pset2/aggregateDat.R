for (dir in 'out') {
  output.files = list.files(dir, full.names=T)
  for(file in output.files) {
    print(file)
  }
}


