#' Dados de tiroteio do Rio de Janeiro
#'
#' @description A site [Fogo Cruzado](https://fogocruzado.org.br/) é uma plataforma digital colaborativa para monitoramento
#' civil da violência armada no Rio de Janeiro e no Recife. Usuários instalam um aplicativo no celular e
#' registram ocorrências de tiroteios (de maneira anônima). Os dados são disponibilizados publicamente através
#' de uma API, feito pela equipe do [`voltdatalab`](https://voltdata.info/). Depois de criar um cadastro no
#' [site](https://api.fogocruzado.org.br/) os dados de tiroteio podem ser extraídos diretamente a API com apoio
#' do pacote `crossfire`. Para tornar os dados no nível município precisar de alguma intervenção manual.
#' Estes dados do pacote estão no nível bairro e já foram tratados e a consulta foi feita em novembro/2020.
#'
#' @name rio_tiroteio
#' @usage data("rio_tiroteio", package = "llsar")
NULL

#' Dados de shape da cidade do Rio de Janeiro
#' @description Objeto da classe `SpatialPolygonsDataFrame`. Os dados de shapes de todos os bairros do Rio de Janeiro (163 bairros, em 2020) para plotar os gráficos, podem ser baixados diteramente do [DataRio](https://www.data.rio/datasets/limite-de-bairros) depois importados pelo R com apoio dos pacotes `rgdal e sp`.
#' @name rio_shape
#' @usage data("rio_shape", package = "llsar")
NULL

#' Matriz simética de proximidades dos bairros cidade do Rio de Janeiro
#' @description Matriz de 0 e 1 com as ajacências dos municípios. Se um município possui adjacência com outro recebe 1 senão, recebe 0.
#'
#' @name rio_proximidades
#' @usage data("rio_proximidades", package = "llsar")
NULL
