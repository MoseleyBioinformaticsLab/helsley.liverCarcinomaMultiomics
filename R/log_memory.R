log_memory = function(){
	linux_memory = system("cat /proc/meminfo", intern = TRUE)
	linux_memory = grep("^MemTotal|^MemAvailable|^Active|^SwapTotal|^SwapFree", linux_memory, value = TRUE)
	linux_memory = grep("anon|active|file", linux_memory, value = TRUE, invert = TRUE)
		
	memory_values = stringr::str_extract(linux_memory, "[[:digit:]].*")
	memory_numbers = as.numeric(stringr::str_extract(memory_values, "[[:digit:]].* "))
	memory_ids = stringr::str_extract(linux_memory, "^[[:alpha:]]+")
	names(memory_numbers) = memory_ids
	memory_string = paste0("Memory: ", paste(paste(c("Total: ", "Available: ", "Active: ", "SwapTotal: ", "SwapFree: "), memory_values, sep = ""), collapse = ", "))
		
	active_to_total = memory_numbers["Active"] / memory_numbers["MemTotal"]
	swapfree_to_swap = memory_numbers["SwapFree"] / memory_numbers["SwapTotal"]
	if (is.nan(swapfree_to_swap)) {
		swapfree_to_swap = 0
	}
		
	swapfree_to_swap = (memory_numbers["SwapTotal"] - memory_numbers["SwapFree"]) / memory_numbers["SwapTotal"]
	if (is.nan(swapfree_to_swap)) {
		swapfree_to_swap = 0
	}
		
		
	if ((active_to_total >= 0.95) || (swapfree_to_swap >= 0.7)) {
		memory_string2 = paste0("HIGH MEMORY USAGE!!! ", memory_string)
		logger::log_warn(memory_string2, namespace = "global")
	} else {
		logger::log_info(memory_string, namespace = "global")
			
	}
}

run_memory_logging = function(log_file, 
															n_iter = Inf,
															interval = 30)
{
	logger::log_appender(logger::appender_tee(log_file), namespace = "global")
	
	i_iter = 0
	while (i_iter <= n_iter) {
		log_memory()
		Sys.sleep(interval)
		i_iter = i_iter + 1
	}
}