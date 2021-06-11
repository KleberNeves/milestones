#' Build the citation network
#'
#' It will call other functions to build the citation
#' net and simplify it into an acyclic graph.
#'
#' @param M A bibliometrix M data frame.
#' @return An igraph citation network.
#' @export
build_network = function (M) {
  cit.net = make_citnet(M)
  cat("\n")
  pajek.net = make_net_for_pajek(cit.net)
  cat("Done!\n")
  pajek.net
}

#' Saves the citation network
#'
#' Saves a citation network in GML and also in a
#' .net format readable by the Pajek software.
#'
#' @param M A bibliometrix M data frame.
#' @param filnema Where to save the files. Extensions will be appended.
#' @return An igraph citation network.
#' @export
save_network <- function (NET, filename) {
  cat("Saving network ...\n")
  igraph::write_graph(NET, paste0(filename,".gml"), format = "GML")
  igraph::write_graph(NET, paste0(filename,".net"), format = "pajek")
}


#' Save list of main path papers
#'
#' Saves the igraph object of a network as a GML file and also
#' saves the list of papers in the main path, as well as the
#' full Web of Science record for the papers in the main path
#' (if given a bibliometric data frame).
#'
#' @param NET The main path network, as an igraph object.
#' @param output_path Path to a directory to save the main path files.
#' @param M A bibliometrix M data frame.
#' @return An bibliometrix data frame with only the main path nodes.
#' @importFrom magrittr %>%
#' @importFrom utils write.table
#' @export
save_papers_mp <- function (NET, filename, M) {
  igraph::write_graph(NET,
                      paste0(filename, "_net.gml"),
                      format = "GML")

  mpPapers = igraph::V(NET)$name
  mpPapers = c(unlist(purrr::map(mpPapers, stringr::str_split, "---")))
  mpPapers = mpPapers[!duplicated(mpPapers)]

  mpM = M[rownames(M) %in% mpPapers,]

  write.table(mpM %>% dplyr::select(SR_FULL, TI, PY, SO, DI, TC),
              paste0(filename, "_papers_list.tsv"),
              sep = "\t", row.names = F)

  write.table(mpM,
              paste0(filename,"_papers_list_full_records.tsv"),
              sep = "\t", row.names = F)
}


#' Pajek file reader: obtains the main path saved from Pajek
#'
#' The main path file will be exported from Pajek. This function
#' reads it back into R as an igraph network built by subsetting
#' the full network.
#'
#' Code is adapted from function *read_net*, from Jonathan H. Morgan (2019)
#' http://mrvar.fdv.uni-lj.si/pajek/R/RMorgan.htm
#'
#' @param fname The path to the Pajek exported file with the main path network.
#' @param net The igraph full network.
#' @return An igraph network with only the main path nodes.
#' @export
read_network <- function (filename, net) {
  read_net <- function(net_file) {
    net <- readLines(net_file)

    trim <- function (x) gsub("^\\s+|\\s+$", "", x)

    #Determining the length of the resulting nodes file
    vertices <- net[[1]]
    vertices <- strsplit(as.character(vertices),' ')
    vertices <- lapply(vertices, function(x) x[x != ""])
    vertices <- as.numeric(vertices[[1]][[2]])

    #Need to account for the fact that meta-information is in the file
    v_num <- vertices + 1

    #Extracting the nodes list from the network file
    nodes <- net[2:v_num]
    nodes <- trim(nodes)

    nodes <- strsplit(as.character(nodes),' ')
    nodes <- lapply(nodes, function(x) x[x != ""])

    shapes <- nodes
    for (i in seq_along(shapes)){
      shapes[[i]] <- nodes[[i]][-c(1:5)]
    }

    shapes <- lapply(shapes, function(x) paste(x,collapse=" "))

    for (i in seq_along(nodes)){
      nodes[[i]] <- nodes[[i]][1:5]
    }

    nodes <- as.data.frame(matrix(unlist(nodes), nrow = length(nodes), byrow = TRUE), stringsAsFactors = FALSE)
    colnames(nodes) <- c('ID', 'Label', 'x-coord', 'y-coord', 'z-coord')

    shapes <-  as.data.frame(matrix(unlist(shapes), nrow = length(shapes), byrow = TRUE),  stringsAsFactors = FALSE)
    colnames(shapes) <- c('shapes information')

    nodes <- cbind(nodes, shapes)

    #Removing shape information if there no information
    nodes[nodes==""] <- NA
    nodes <- nodes[, colSums(is.na(nodes)) == 0]

    rm(shapes, i)

    #Creating Edges File
    vertices <- (vertices + 1)

    edges <- net[-c(1:vertices)]

    `Tie Type` <- edges[[1]]

    edges <- edges[-c(1)]
    edges <- trim(edges)

    edges <- strsplit(as.character(edges),' ')
    edges <- lapply(edges, function(x) x[x != ""])

    edges <-  as.data.frame(matrix(unlist(edges), nrow = length(edges), byrow = TRUE), stringsAsFactors = FALSE)
    colnames(edges) <- c('Person i', 'Person j', 'Weight')

    #Removing edge color information as it adds little value when importing data
    edges <- edges[-c(4:5)]

    edges$`Tie Type` <- `Tie Type`

    rm(net, `Tie Type`)

    #writing objects to the Global Environment
    return (list(nodes = nodes, edges = edges))
  }

  x = read_net(filename)
  vertices = x$nodes
  ties = x$edges
  rm(x)

  vertices$Node = as.numeric(gsub('v|\"', "", vertices$Label))

  ties$Source = vertices$Node[as.numeric(ties$`Person i`)]
  ties$Target = vertices$Node[as.numeric(ties$`Person j`)]
  main_ties = c(rbind(ties$Source,ties$Target))

  main_path_nodes = igraph::V(net)[vertices$Node]
  main_path_edges = igraph::get.edge.ids(net, vp = main_ties, directed = T)

  main_path = igraph::subgraph.edges(net, main_path_edges, delete.vertices = T)
  main_path
}

#' Matches citations to build citation net
#'
#' Matches the names in the cited references field to the
#' short titles to identify papers in the dataset that cited
#' each other. It's a version of the original bibliometrix
#' version with a few steps omitted for performance.
#'
#' @param M The bibliometrix data frame.
#' @param min.citations The minimum number of received citations for the paper to be included.
#' @return A list with the citation information and metadata.
make_hist_citation_net <- function (M, min.citations = 1) {
  min.citations = max(c(1, min.citations))
  M$TC = as.numeric(M$TC)
  M = M[!is.na(M$TC), ]
  if (!("SR_FULL" %in% names(M))) {
    M = bibliometrix::metaTagExtraction(M, Field = "SR")
  }
  M = M[order(M$PY), ]
  M2 = M[M$TC >= min.citations, ]
  if (dim(M2)[1] == 0) {
    cat("\nNo document has a number of citations above the fixed threshold\n")
    return(NULL)
  }
  N = dim(M2)[1]
  N2 = dim(M)[1]
  rows = c(1:N2)
  lCit = Matrix::Matrix(0, N, N2)

  for (i in 1:N) {
    if (i %% 100 == 0 | i == N) cat(paste0("Articles processed: ", i, "\n"))

    # Searches by name
    x = M2$SR_FULL[i]
    Year = M2$PY[i]
    pos = stringr::str_which(M$CR[M$PY >= Year], x)
    pos = rows[M$PY >= Year][pos]

    # Searches by DOI
    # if (!is.na(M2$DI[i])) {
    #   pos2 = stringr::str_which(fixed(M$CR[M$PY >= Year]), M2$DI[i])
    #   pos2 = rows[M$PY >= Year][pos2]
    #   p1 = pos
    #   pos = unique(pos, pos2)
    #   if (length(p1) > 0 & length(pos) > 0) {
    #     if (p1 != pos) { cat("AAAAAAAAA")}
    #   }
    # }

    if (length(pos) > 0) {
      lCit[i, pos] = 1
    }
  }

  LCS = Matrix::rowSums(lCit)
  ind = which(LCS > M2$TC)
  LCS[ind] = M2$TC[ind]
  M2$LCS = LCS
  row.names(lCit) = M2$SR
  colnames(lCit) = M$SR
  lCit = lCit[, (M$SR %in% M2$SR)]
  if (!("DI" %in% names(M2))) {
    M2$DI = NA
  }
  df = data.frame(Paper = M2$SR, DOI = M2$DI, Year = M2$PY,
                  LCS = LCS, GCS = M2$TC, stringsAsFactors = F)
  df = df[order(df$Year), ]
  row.names(df) = paste(df$Year, rep("-", dim(df)[1]), 1:dim(df)[1])
  results = list(NetMatrix = Matrix::t(lCit), histData = df, M = M2,
                 LCS = LCS)
  return(results)
}

#' Builds a historical citation igraph net
#'
#' Uses a M bibliometrix data frame to build the citation network
#' for the papers in the dataset, in the igraph format.
#'
#' @param M The bibliometrix data frame.
#' @return A igraph network made of citations between papers.
#' @export
make_citnet <- function(M) {
  cat("Building edges ...\n")
  histResults = make_hist_citation_net(M, min.citations = 1)
  cat("Making graph ...\n")
  ADJ = as.matrix(histResults$NetMatrix)
  NET = igraph::graph_from_adjacency_matrix(ADJ, mode = "directed", diag = F)
  cat("\nSimplifying network ...")
  NET = igraph::simplify(NET, remove.multiple = T, remove.loops = T)
  NET
}

#' Removes loops from a network so that Pajek can read it
#'
#' Removes loops of order up to 4 to simplify the network
#' and make it a direct acyclic graph (DAG). This is necessary
#' so that Pajek can run the main path analysis algorithm.
#'
#' Called for the side effects: it will save GML and Pajek
#' files of the network (Historical Citation Net). It will
#' warn you if the resulting network is still not a DAG.
#'
#' @param NET An igraph network.
#' @return The simplified network.
#' @export
make_net_for_pajek <- function (NET) {
  NET = igraph::simplify(NET, remove.multiple = T, remove.loops = T)
  # browser()
  tries = 0
  while (!igraph::is_dag(NET) & tries < 10) {
    removed = F
    for (size in 2:4) {
      w = rep(0, 4)
      m = igraph::as_adjacency_matrix(NET, sparse = T)
      mm = m
      for (xp in 2:4) {
        mm = mm %*% m
        w[xp] = sum(Matrix::diag(mm))
      }

      if (w[size] > 0) {
        cat(paste0("\nFinding citation loops of size ", size," ...\n"))
        loops = find_cycles(NET, size)

        cat(paste0("Found ", length(loops), " loops."))
        if (length(loops) > 0) {
          cat(paste0("\nRemoving citation loops of size ", size," ...\n"))
          NET = make_loops_into_families(loops, NET)

          NET = igraph::simplify(NET, remove.multiple = T, remove.loops = T)
          NET = igraph::delete.vertices(NET, igraph::degree(NET) == 0)
          removed = T
          break
        }
      }
    }
    if (!removed) {
      tries = tries + 1
    }
  }

  if (!igraph::is_dag(NET)) {
    warning("Warning: network still contains loops (not a DAG).")
  }

  NET
}

#' Find cycles of a specified length in a graph
#'
#' Helper function for making networks acyclical. It
#' finds loops of a given order in an igraph network.
#'
#' @param g The igraph network.
#' @param size The length of the loops to find.
#' @return A list of vectors, each vector providing the ids of the nodes in the loop.
find_cycles <- function(g, size) {
  Cycles = NULL
  if (size == 2) {
    for(v1 in igraph::V(g)) {
      if (v1 %% 100 == 0) cat(paste0("Articles processed: ", v1, "\n"))
      nei1 = igraph::neighbors(g, v1, mode="out")
      nei1 = nei1[nei1 > v1]
      for(v2 in nei1) {
        nei2 = igraph::neighbors(g, v2, mode="out")
        if (v1 %in% nei2) {
          Cycles = c(Cycles, list(c(v1,v2)))
        }
      }
    }
    Cycles = unique(
      lapply(Cycles, function (x) {
        unique(sort(unlist(x)))
      })
    )
  } else if (size == 3) {
    for(v1 in igraph::V(g)) {
      if (v1 %% 100 == 0) cat(paste0("Articles analyzed ", v1, "\n"))
      nei1 = igraph::neighbors(g, v1, mode="out")
      nei1 = nei1[nei1 > v1]
      for(v2 in nei1) {
        nei2 = igraph::neighbors(g, v2, mode="out")
        nei2 = nei2[nei2 > v2]
        for(v3 in nei2) {
          nei3 = igraph::neighbors(g, v3, mode="out")
          if (v1 %in% nei3) {
            Cycles = c(Cycles, list(c(v1,v2,v3)))
          }
        }
      }
    }
    Cycles = unique(
      lapply(Cycles, function (x) {
        unique(sort(unlist(x)))
      })
    )
  } else {
    for(v1 in igraph::V(g)) {
      if (v1 %% 100 == 0) cat(paste0("Articles analyzed ", v1, "\n"))
      nei1 = igraph::neighbors(g, v1, mode="out")
      nei1 = nei1[nei1 > v1]
      for(v2 in nei1) {
        nei2 = igraph::neighbors(g, v2, mode="out")
        nei2 = nei2[nei2 > v2]
        for(v3 in nei2) {
          nei3 = igraph::neighbors(g, v3, mode="out")
          nei3 = nei3[nei3 > v3]
          for(v4 in nei3) {
            nei4 = igraph::neighbors(g, v4, mode="out")
            if (v1 %in% nei4) {
              Cycles = c(Cycles, list(c(v1,v2,v3,v4)))
            }
          }
        }
      }
    }
    Cycles = unique(
      lapply(Cycles, function (x) {
        unique(sort(unlist(x)))
      })
    )
  }
  Cycles
}

#' Combines attribute for merged paper nodes
#'
#' When simplifying networks by merging nodes/papers, this
#' is the function used to combine the paper titles in one.
#'
#' @param x A vector with the titles of the papers.
#' @return The combined titles.
family_attr_comb <- function(x) {
  paste(x, collapse = "---")
}

#' Merges papers which are in a citation loop
#'
#' Merges the papers passed to it into one node. This is used
#' when removing loops to make the citation net acyclic. If A
#' cites B and B cites A, this will become node AB, with their
#' outgoing and ingoing edges pooled together.
#'
#' @param loops A list of loops, each a vector with node ids.
#' @param g The citation igraph network.
#' @return A data frame with the information extracted.
make_loops_into_families <- function (loops, g) {
  contraction = 1:igraph::vcount(g)
  contracted = numeric(0)

  for (loop in loops) {
    contraction[loop] = loop[1]
    contracted = c(contracted, loop[2:length(loop)])
  }
  g = igraph::contract(g, mapping = contraction, vertex.attr.comb = family_attr_comb)
  contracted = contracted[contracted %in% igraph::V(g)]
  g = igraph::delete.vertices(g, igraph::V(g)$name[contracted])

  g
}

