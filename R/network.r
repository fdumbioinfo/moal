#' @title stringDB interaction network
#' @param nodelist character list of Symbols to make network
#' @param foldchange data.frame Symbol list with fold-change values to display
#' @param pval data.frame Symbol list with  p-values to display
#' @param species character for species hs mm rn dr
#' @param layout numeric for layout neetwork 1 fr by default 2 dh 3 tree 4 circle 5 grid 6 sphere
#' @param seed numeric to change random layout
#' @param intmaxdh numeric maximum number of interaction to use for Davidson and Harel algorithm layout
#' @param nodelabelsize numeric change node label size
#' @param nodesize numeric change edge size
#' @param edgeweight numeric change node edge weight
#' @param edgewidth numeric change node edge weight
#' @param title character plot title
#' @param path character file path
#' @param filename character results pathways file name
#' @return result folder
#' @examples
#' # not run
#' # network(symbollist)
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by left_join select
#' @importFrom rlang .data
#' @importFrom stats setNames
#' @importFrom igraph graph_from_data_frame degree V layout_with_fr layout_with_dh layout_as_tree layout_in_circle layout_on_grid layout_on_sphere 
#' @importFrom foreach foreach
#' @importFrom colourvalues colour_values
#' @importFrom grDevices colorRamp
#' @importFrom Rgraphviz plot
#' @import moalstringdbhs
#' @import moalstringdbmm
#' @import moalstringdbrn
#' @import moalstringdbdr
#' @import moalstringdbss
#' @noRd
network <- function(
    nodelist = NULL, foldchange = NULL, pval = NULL, species = NULL, layout = 1, seed = 1234567,
    intmaxdh = 5000, nodelabelsize = 0.39, nodesize = 0.5, edgeweight = 0.1, edgewidth = 0.2, title = "Network", path = ".", filename = NULL)
{
  i=j=k=m=1
  c("layout_with_fr","layout_with_dh","layout_as_tree","
    layout_in_circle","layout_on_grid","layout_on_sphere") -> Layout0
  Layout0[layout] -> Layout1
  # species selection
  ifelse(
    is.null(species),
    orthoinfo[[6]] -> Species0,
    orthoinfo %>% sapply("[[",1) %>% grep(species,.) %>% "[["(orthoinfo,.) -> Species0)
  # Nodes
  nodelist %>% as.character -> NodeList0
  nodelist %>% as.character -> NodeList1
  NodeList0 %>% annot( species=Species0[1], dboutput="ebi", idtype = "SYMBOL" ) -> NodeList1
  # check
  if( nrow(NodeList1) > 0 &
      (((((NodeList1$ENSGID %>% is.na)|(NodeList1$ENSGID %>% "=="(.,""))|(NodeList1$ENSGID %>% "=="(.,"NA"))) %>% all)|
      (((NodeList1$ENSPIDs %>% is.na)|(NodeList1$ENSPIDs %>% "=="(.,""))|(NodeList1$ENSPIDs %>% "=="(.,"NA"))) %>% all)) %>%
      any %>% "!"(.)))
    { NWNODE1 <- TRUE }else{ NWNODE1 <- FALSE }
  #
  if(NWNODE1)
  {
    # Nodes
    NodeList1$ENSPIDs %>% strsplit("\\|") %>% unlist %>% unique -> NodeList2
    NodeList2 %>% is.na %>% which -> sel ; if(length(sel) > 0){ NodeList2[-sel] -> NodeList2 }
    NodeList2 %>% "=="(.,"NA") %>% which -> sel ; if(length(sel) > 0){ NodeList2[-sel] -> NodeList2 }
    NodeList2 %>% as.character %>% moal::annot(species=Species0[1], idtype="ENSP") %>%
      dplyr::select(.data$ENSPID,.data$Symbol) %>% data.frame -> NodeList3
    # check
    if(nrow(NodeList3) > 0 & NodeList3$Symbol %>% is.na %>% all %>% "!"(.)){ NWNODE3 <- TRUE }else{ NWNODE3 <- FALSE }
    #
    if(NWNODE3)
    {
      # stringdb loading
      if( Species0[1] == "hs" ){ moalstringdbhs::stringdbhs -> EdgeList0 }
      if( Species0[1] == "mm" ){ moalstringdbmm::stringdbmm -> EdgeList0 }
      if( Species0[1] == "rn" ){ moalstringdbrn::stringdbrn -> EdgeList0 }
      if( Species0[1] == "dr" ){ moalstringdbdr::stringdbdr -> EdgeList0 }
      if( Species0[1] == "ss" ){ moalstringdbss::stringdbss -> EdgeList0 }
      # Edges
      NodeList3 %>% dplyr::inner_join(EdgeList0, by =c("ENSPID"="NodeA")) %>% setNames(c("NodeA","SymbolA","NodeB","CombinedScore")) -> EdgeList1
      # Direct connection filtering
      NodeList3 %>% dplyr::inner_join(EdgeList1, by=c("ENSPID"="NodeB")) %>% dplyr::select(c(1:3,5,4)) %>% setNames(c("NodeA","SymbolA","NodeB","CombinedScore","SymbolB")) -> EdgeList2
      # Confidence filtering
      seq(0,999,10) -> Confidence0
      mode(EdgeList2$CombinedScore) <- "numeric"
      EdgeConfNumber <- foreach(k=1:length(Confidence0)) %do%
        { EdgeList2 %>% dplyr::filter( .data$CombinedScore >= Confidence0[k] ) %>% nrow }
      #
      if(length(EdgeConfNumber %>% unlist %>% "<"(intmaxdh) %>% which) > 0)
      {
        EdgeConfNumber %>% unlist %>% "<"(intmaxdh) %>% which %>% min %>% Confidence0[.] -> Confidence1
        EdgeList2 %>% dplyr::filter( .data$CombinedScore >= Confidence1 ) -> EdgeList3
        # remove NA
        is.na(EdgeList3$SymbolA) %>% which -> sel ; if( length(sel) > 0 ){ EdgeList3[-sel,] -> EdgeList3  }
        EdgeList3$SymbolA %>% grep("^$|^NA$|^row",.) -> sel ; if( length(sel) > 0 ){ NodeList3[-sel,] -> NodeList3  }
        is.na(EdgeList3$SymbolB) %>% which -> sel ; if( length(sel) > 0 ){ EdgeList3[-sel,] -> EdgeList3  }
        EdgeList3$SymbolB %>% grep("^$|^NA$|^row",.) -> sel ; if( length(sel) > 0 ){ NodeList3[-sel,] -> NodeList3  }
        if( nrow(EdgeList3) > 0 & nrow(EdgeList3) < intmaxdh )
        {
          EdgeList3 %>% dplyr::select(.data$NodeA,.data$NodeB) -> EdgeList4
          EdgeList3 %>% dplyr::select(.data$NodeA,.data$SymbolA) %>% setNames(c("ENSPID","Symbol")) %>% data.frame ->  NodeListA
          EdgeList3 %>% dplyr::select(.data$NodeB,.data$SymbolB) %>% setNames(c("ENSPID","Symbol")) %>% data.frame ->  NodeListB
          rbind(NodeListA,NodeListB) %>% unique -> NodeList4
          #
          igraph::graph_from_data_frame(d=EdgeList4, vertices=NodeList4, directed=F) -> NetWork0
          igraph::degree(NetWork0, mode="all") -> Deg0
          # layout
          # NetworkLayout1 <- "layout_with_dh"
          # NetworkLayout1 <- layout
          NetworkLayout1 <- Layout1
          # Nodes
          # data.frame(Symbol=V(NetWork0)$Symbol %>% as.character,ENSPID=V(NetWork0) %>% attr("names") %>% as.character,stringsAsFactors=F) %>%
          #   inner_join(NodeList1) %>% data.frame -> NodeList
          # Edges output
          if(is.null(filename)){ filename <- "network" }
          paste("Edges_",filename,"_",nrow(EdgeList3),"_",length(V(NetWork0)$Symbol),".tsv",sep ="") -> FileName0
          EdgeList3 %>% output(file.path(path,FileName0))
          #
          set.seed(seed)
          eval(parse(text=paste(NetworkLayout1[[i]][1],"(NetWork0)",sep =""))) -> LayOut
          #
          rep("gray",V(NetWork0)$Symbol %>% length) -> ColorNode
          rep("gray",V(NetWork0)$Symbol %>% length) -> FrameColorNode
          ### fold-change
          if(!is.null(foldchange))
          {
            # color range up
            grDevices::colorRampPalette(c("green4","green","greenyellow")) -> ColorRangeDown0
            ColorRangeDown0(5) -> ColorRangeDown1
            ColorRangeDown1 %>% data.frame(color=.,value=c(-5:-1)) -> ColorRangeDown2
            # color range down
            grDevices::colorRampPalette(c("yellow","orange","red")) -> ColorRangeUp0
            ColorRangeUp0(5) -> ColorRangeUp1
            ColorRangeUp1 %>% data.frame(color=.,value=c(1:5)) -> ColorRangeUp2
            # rbind(ColorRangeDown2,ColorRangeUp2) -> ColorRange0
            # ColorRange0
            #
            foreach(j=1:length(foldchange)) %do%
              {
                V(NetWork0)$Symbol %>% data.frame(Symbol=.) %>% dplyr::left_join(foldchange[[j]]) -> FoldChange0
                # foldchange[[j]] -> tt
                if( nrow(FoldChange0) > 1 )
                {
                  FoldChange0 %>% dplyr::arrange(.data[[colnames(FoldChange0)[2]]]) -> FoldChange1
                  FoldChange0 -> FoldChange1
                  rep("gray",nrow(FoldChange1)) -> ColorNode0
                  # FoldChange1 %>% dplyr::select(.data[[colnames(FoldChange1)[2]]]) %>% unlist %>% "<"(ColorRange0$value[1]) %>% which -> selc1
                  FoldChange1 %>% dplyr::select(.data[[colnames(FoldChange1)[2]]]) %>% unlist %>% "<"(ColorRangeDown2$value[1]) %>% which -> selc1
                  if(length(selc1) > 0){ ColorNode0 %>% replace(selc1,ColorRangeDown2$color[1]) -> ColorNode0 }
                  # selc19 <- foreach(m=1:(nrow(ColorRange0)-1)) %do%
                  #   {
                  #     ((FoldChange1 %>% dplyr::select(.data[[colnames(FoldChange1)[2]]]) %>% unlist %>% ">"(ColorRange0$value[m+1])) &
                  #        (FoldChange1 %>% dplyr::select(.data[[colnames(FoldChange1)[2]]]) %>% unlist %>% "<"(ColorRange0$value[m+1+1]))) %>% which
                  #   }
                  # foreach(m=1:length(selc19)) %do% { ColorNode0 %>% replace(selc19[[m]],ColorRange0$color[m]) -> ColorNode0 }
                  # ColorRangeDown2
                  selc19down <- foreach(m=1:(nrow(ColorRangeDown2)-1)) %do%
                    {
                      ((FoldChange1 %>% dplyr::select(.data[[colnames(FoldChange1)[2]]]) %>% unlist %>% ">="(ColorRangeDown2$value[m])) &
                         (FoldChange1 %>% dplyr::select(.data[[colnames(FoldChange1)[2]]]) %>% unlist %>% "<"(ColorRangeDown2$value[m+1]))) %>% which
                    }
                  foreach(m=1:length(selc19down)) %do% { ColorNode0 %>% replace(selc19down[[m]],ColorRangeDown2$color[m+1]) -> ColorNode0 }
                  # ColorRangeUp2
                  selc19up <- foreach(m=1:(nrow(ColorRangeUp2)-1)) %do%
                    {
                      ((FoldChange1 %>% dplyr::select(.data[[colnames(FoldChange1)[2]]]) %>% unlist %>% ">="(ColorRangeUp2$value[m])) &
                         (FoldChange1 %>% dplyr::select(.data[[colnames(FoldChange1)[2]]]) %>% unlist %>% "<"(ColorRangeUp2$value[m+1]))) %>% which
                    }
                  foreach(m=1:length(selc19up)) %do% { ColorNode0 %>% replace(selc19up[[m]],ColorRangeUp2$color[m+1]) -> ColorNode0 }
                  #
                  FoldChange1 %>% dplyr::select(.data[[colnames(FoldChange1)[2]]]) %>% unlist %>% ">"(ColorRangeUp2$value[5]) %>% which -> selc10
                  if(length(selc10) > 0){ ColorNode0 %>% replace(selc10,ColorRangeUp2$color[5]) -> ColorNode0 }
                  FoldChange1$Symbol %>% data.frame(Symbol=.,ColorNode=ColorNode0,check.names = F) -> FoldChange2
                  # V(NetWork0)$Symbol %>% data.frame(Symbol=.) %>% right_join(FoldChange2) %>% dplyr::select(.data$ColorNode) %>% unlist -> ColorNode
                  FoldChange2$ColorNode -> ColorNode
                  ColorNode -> FrameColorNode
                }
                if(!is.null(pval))
                {
                  V(NetWork0)$Symbol %>% data.frame(Symbol=.) %>% dplyr::left_join(pval[[j]]) -> Pval0
                  c(0.05,0.01,0.001) -> tr0
                  c("*","**","***") -> trs0
                  Pval0$Symbol -> PvalNode
                  foreach(i=1:length(tr0)) %do%
                    {
                      Pval0 %>% dplyr::select(.data[[colnames(Pval0)[2]]]) %>% unlist %>% "<"(tr0[i]) %>% which -> sel
                      PvalNode[sel] %>% paste(trs0[i],sep="") -> values
                      PvalNode[sel] %>% paste("*",sep="") -> values
                      PvalNode %>% replace(sel,values) -> PvalNode
                      FrameColorNode %>% replace(sel,"purple") -> FrameColorNode 
                    }
                  # 
                }
                #
                ifelse(!is.null(foldchange),colnames(FoldChange0)[2] %>% gsub("fc_","",.) -> compName, "nofc" -> compName)
                ColorNode %>% grep("gray",.,invert = T) %>% length -> FCNode0
                ColorNode %>% length -> AllNode0
                paste(filename,"_",compName,"_",nrow(EdgeList3),"_",FCNode0,"_",AllNode0,".pdf",sep="") -> FileName0
                FileName0 %>% gsub("(.*).pdf","\\1",.) -> Title0
                file.path(path,FileName0) -> FileName1
                #
                pdf(FileName1 , width = 10 , height = 10 )
                # plot( NetWork0, rescale=T, main=Title0, cex.main=1,
                Rgraphviz::plot( NetWork0, rescale=T,
                      edge.width=edgewidth, edge.weight=edgeweight, edge.color="gray",
                      layout=LayOut, edge.curved=0.1, arrow.mode=3, edge.arrow.size=0.1, edge.arrow.width=0.5,
                      vertex.label=PvalNode, vertex.size=(log2(Deg0)+5)*nodesize, vertex.label.cex=nodelabelsize,
                      vertex.frame.color=FrameColorNode, vertex.label.font=2, vertex.label.color="blue", vertex.color=ColorNode )
                title(main=Title0,cex.main=0.7)
                # legend(1, 95, legend=c("Line 1", "Line 2"),
                #        col=c("red", "blue"), lty=1:2, cex=0.8)
                graphics.off()
              }
          }else{
            # ifelse(!is.null(foldchange),colnames(foldchange)[2] %>% gsub("fc_","",.) -> compName, "nofc" -> compName )
            ColorNode %>% grep("gray",.,invert = T) %>% length -> FCNode0
            ColorNode %>% length -> AllNode0
            # paste(filename,"_",compName,"_",nrow(EdgeList3),"_",FCNode0,"_",AllNode0,".pdf",sep="") -> FileName0
            paste(filename,"_",nrow(EdgeList3),"_",FCNode0,"_",AllNode0,".pdf",sep="") -> FileName0
            FileName0 %>% gsub("(.*).pdf","\\1",.) -> Title0
            file.path(path,FileName0) -> FileName1
            #
            pdf(FileName1 , width = 10 , height = 10 )
            Rgraphviz::plot( NetWork0, rescale=T, main=Title0, cex.main=2,
                  edge.width=edgewidth, edge.weight=edgeweight, edge.color="gray",
                  layout=LayOut, edge.curved=0.1, arrow.mode=3, edge.arrow.size=0.1, edge.arrow.width=0.5,
                  vertex.label=V(NetWork0)$Symbol, vertex.size=(log2(Deg0)+5)*nodesize, vertex.label.cex=nodelabelsize,
                  vertex.frame.color=FrameColorNode, vertex.label.font=2, vertex.label.color="blue", vertex.color=ColorNode )
            legend(1, 95, legend=c("Line 1", "Line 2"),
                   col=c("red", "blue"), lty=1:2, cex=0.8)
            graphics.off()
          }
        }
      }
    }# fi NWNODE3
  }# fi NWNODE1
}
