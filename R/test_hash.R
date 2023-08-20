#' @useDynLib perms test_hash


#' @export
Rtest_hash= function(){
  
  res = .Call(test_hash)
  
  return
}
