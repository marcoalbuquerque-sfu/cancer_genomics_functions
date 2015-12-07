read_biomart_gene_file <- function(biomart_gene_file) {

    data <- read.table(
        file=biomart_gene_file,
        header=F,
       	stringsAsFactors=F
    	)

    data <- data[,c(5,4,3,6)]

    colnames(data) <- c(
        "chr",
        "start",
        "end",
        "gene"
    	)

    data <- data[order(data$end),]
    data <- data[order(data$start),]
    data <- data[order(data$chr),]
    
    return(data)

    }

read_titan_cnv_file <- function(titan_cnv_file) {

    data <- read.table(
        file=titan_cnv_file,
        header=T,
        stringsAsFactors=F
    	)

    data <- data[,c(2,3,4,1,7,9,10)]

    colnames(data) <- c(
        "chr",
        "start",
        "end",
        "sample",
        "logR",
        "state",
        "copy"
    	)

    return(data)

    }

generate_common_genome_frame <- function() {

    data <- data.frame(
        chr = c(1:22, "X", "Y"),
        end = c(
			249250621,
			243199373,
			198022430,
			191154276,
			180915260,
			171115067,
			159138663,
			146364022,
			141213431,
			135534747,
			135006516,
			133851895,
			115169878,
			107349540,
			102531392,
			90354753,
			81195210,
			78077248,
			59128983,
			63025520,
			48129895,
			51304566,
			155270560,
			59373566
			)
        )

    return(data)

    }

calculate_percent_genome_alteration <- function(cnv_frame, genome_frame) {
    
    truth <- !cnv_frame$copy == 2;
    numer <- sum( as.numeric( cnv_frame[truth,]$end - cnv_frame[truth,]$start) );
    denom <- sum( as.numeric( genome_frame$end - 1 ) );
    data <- data.frame(
        sample = cnv_frame$sample[1],
    	pga = numer / denom,
    	stringsAsFactors=F
    	)
    return(data) 

    }

read_mutsigCV_genes <- function(mutsigCV_file) {

    data <- read.table(
        mutsigCV_file,
        stringsAsFactors=F,
        header=T
        )

    data <- data[,c(1,14,15)]

    colnames(data) <- c(
        "gene",
        "p.value",
        "q.value"
        )

    return(data)

    }

read_oncodriveFM_genes <- function(oncodriveFM_file) {

    data <- read.table(
    	oncodriveFM_file,
    	stringsAsFactors=F,
        header=T
    	)

    data <- data[,c(1,2,3)]

    data <- data.frame(
        "gene"=rownames(data),
        "p.value"=data[,1],
        "q.value"=data[,2],
        stringsAsFactors=F
    	)

    return(data)

    }

filter_gene_list <- function(gene_list, q.value=0.05, p.value=1) {
    
    data <- gene_list;
    data <- data[data$q.value <= q.value,];
    data <- data[data$p.value <= p.value,];
    return(data)

    }

join_gene_list <- function(mutsig=NULL, oncodrive=NULL, other=NULL) {

	data <- NULL

	if (!is.null(mutsig)) {
        data <- mutsig;
	    }

	if (!is.null(oncodrive)) {
        if (is.null(data)) {
            data <- oncodrive;
            }
        else {
            data <- rbind(data, oncodrive);
            }
	    }

	if (!is.null(other)) {
        if (is.null(data)) {
            data <- other;
            }
        else {
            data <- rbind(data, other);
            }
	    }

	data <- data[!duplicated(data$gene),]

	return(data)

    }

read_titan_params_file <- function(params_file) {

    data <- read.table(
        params_file,
        fill = TRUE,
        stringsAsFactors=FALSE
        )
   
    return(data)

    }

fetch_validity_index <- function(params_frame) {
    return(as.numeric(params_frame[17,5]))
    }