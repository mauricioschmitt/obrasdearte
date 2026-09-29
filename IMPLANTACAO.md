# Como colocar no ar

O portal são duas páginas, e nenhuma delas precisa de servidor:

- `docs/index.html` — o portal. Clica no mapa, calcula a vazão.
- `docs/preparar.html` — onde você arrasta os GeoTIFFs e ela gera os dados
  de uma região nova. Roda inteira no navegador, sem Python.

```
docs/
  index.html           o portal
  preparar.html        a página de preparo
  dados/               o que a página de preparo gera
hidro/                 versão em Python do mesmo motor (opcional)
preparar_regiao.py     versão em linha de comando (opcional)
```

O `index.html` funciona sozinho. Se não encontrar `dados/manifesto.json`,
abre em modo demonstração com relevo sintético. Dá para publicar hoje e ir
acrescentando regiões reais depois.

---

## Acrescentar uma região

Abra `preparar.html`, no endereço publicado ou dando dois cliques no arquivo
na sua máquina. Funciona dos dois jeitos, e não precisa de QGIS nem de Python.

1. Arraste o GeoTIFF do modelo de elevação. Pode estar em graus: a página
   reprojeta sozinha para UTM, detectando a zona pela longitude.
2. Arraste sobre a prévia do relevo para recortar a região, e escolha a
   resolução. A página mostra a área, o número de células e o tamanho
   estimado do download, para você decidir antes de processar.
3. Arraste o GeoTIFF do MapBiomas. Pode estar em outra resolução, outro
   recorte e outro CRS.
4. Preencha o identificador, o nome e os coeficientes da curva IDF.
5. Clique em "Preparar região" e leia o registro que aparece.
6. Baixe os três arquivos e suba no GitHub, em `docs/dados/`, por
   "Add file" → "Upload files".

A página busca sozinha o `manifesto.json` já publicado e preserva as regiões
que você subiu antes, então dá para ir acrescentando uma por vez.

### Onde baixar os dados sem QGIS

**Modelo de elevação.** O OpenTopography tem uma interface web onde você
desenha o retângulo no mapa, escolhe o Copernicus GLO-30 e baixa um GeoTIFF.
Vem em graus, e tudo bem: a página de preparo reprojeta. O FABDEM tende a
ser melhor em SC porque remove vegetação, e a Mata Atlântica introduz vários
metros de erro em MDE de satélite bruto; confira a licença.

**Uso e cobertura.** MapBiomas, pela área de downloads do site deles, que
oferece recortes por estado e por município.

Não tenho acesso a busca, então confirme endereços, licenças e formatos
antes de fechar o plano de dados.

### Precisão da reprojeção

A reprojeção feita no navegador foi comparada com a do GDAL na mesma grade,
num terreno com 600 m de amplitude e cerca de 10 m de desnível entre células
vizinhas. Erro médio de 0,025 m contra 0,004 m do GDAL, erro absoluto médio
de 3,73 m contra 3,80 m. São equivalentes, e o resíduo vem da dupla
reamostragem, não do método.

A direção de fluxo e as classes de uso do solo batem em 100% com a versão
Python, célula a célula.

### O script em Python continua existindo

`preparar_regiao.py` faz exatamente a mesma coisa pela linha de comando.
A saída dos dois foi comparada célula a célula: a direção de fluxo e as
classes de uso do solo batem em 100%. Use o script quando for processar
muitas regiões de uma vez; para uma ou duas, a página é mais simples.

---

## Processando em lote pela linha de comando

```bash
pip install numpy scipy rasterio pyproj
# opcional, deixa o preenchimento de depressões ~20× mais rápido:
pip install richdem

python preparar_regiao.py \
    --mde recortes/itajai_mirim_mde.tif \
    --uso recortes/itajai_mirim_mapbiomas.tif \
    --id itajai-mirim \
    --nome "Bacia do Itajaí-Mirim" \
    --idf 1132,0.159,12.0,0.771 \
    --saida docs/dados
```

Gera três arquivos em `docs/dados/`: o `.bin.gz` que o usuário baixa, um
`.json` com os metadados, e o `manifesto.json` atualizado. Rodar de novo com
outro `--id` acrescenta região; com o mesmo `--id`, substitui.

### Recorte: por bacia, nunca por município

É a decisão mais importante do projeto. Uma bacia hidrográfica nunca
atravessa um divisor de águas, então recortando por unidade hidrográfica
toda microbacia cai inteira dentro de um arquivo. Recortando por município,
as bacias que cruzam a divisa ficam cortadas ao meio e a área sai menor que
a real, sem nenhum aviso. É um erro silencioso, que é o pior tipo.

Use as ottobacias da Base Hidrográfica Ottocodificada da ANA, ou as bacias
hidrográficas oficiais de SC, subdivididas até chegar no tamanho certo.

**Tamanho alvo: 1.500 a 2.500 km² por arquivo.** A conta é de 1,4 byte por
célula depois da compressão. A 30 m isso dá de 1,5 a 3 MB por download, que
é o teto razoável para a conexão de uma prefeitura pequena. As 10 regiões
hidrográficas de SC inteiras são grandes demais (uma RH daria uns 15 MB);
cada uma precisa virar 3 a 5 arquivos.

### De onde vêm os dados

**Modelo de elevação.** Comece por um MDE de 30 m de cobertura nacional, e
não pelo MDT métrico do SIGSC. Motivo prático: o mesmo script passa a
funcionar em qualquer lugar do Brasil sem mudar nada, que é o objetivo de
vocês. Candidatos: Copernicus GLO-30 e FABDEM. O FABDEM tende a ser melhor
em SC porque remove vegetação e edificação, e a Mata Atlântica introduz
vários metros de erro em MDE de satélite bruto. Confira a licença antes,
porque varia. O MDT do SIGSC vale para um piloto urbano, onde 30 m é grosso
demais, mas gera arquivos pesados demais para uso estadual.

**Uso e cobertura.** MapBiomas, recortado na mesma área. Não precisa
reprojetar nem reamostrar antes: o script alinha sozinho sobre a grade do
MDE, por vizinho mais próximo, que é o correto para dado temático.

Confira os códigos de classe contra a legenda da coleção que você baixar.
Eles mudam entre coleções. O script lista os códigos que não soube traduzir
e quanto de área eles ocupam. Enquanto aparecer código relevante nessa
lista, o CN está errado.

**Curva IDF.** O `--idf` recebe K, a, b, c da equação
`i = K · TR^a / (t + b)^c`, com i em mm/h e t em minutos. Para SC existem
equações de chuvas intensas municipalizadas publicadas pela Epagri; procure
o trabalho de desagregação de chuvas do estado. Os coeficientes do exemplo
acima são inventados, só para o teste. **Não use em projeto.**

Não tenho acesso a busca, então confirme disponibilidade, licença e formato
dessas bases antes de fechar o plano de dados.

---

## Publicar

Copie a pasta `docs/` inteira para qualquer hospedagem de arquivo estático.
Não há servidor, banco nem build.

**Cloudflare Pages** é a recomendação. Crie a conta, escolha "Upload
assets", arraste a pasta `docs/`. Sai um endereço em minutos. Banda
ilimitada, gratuito, sem cartão. Para atualizar, arraste de novo.

**Netlify Drop** (`app.netlify.com/drop`) é ainda mais direto: arrasta a
pasta sem nem criar conta antes. Bom para mostrar em reunião. Migre para
Cloudflare quando virar oficial.

Depois vale pedir um subdomínio à TI da universidade e apontar para lá. O
endereço institucional dá credibilidade com prefeitura e não custa nada.

### Por que não um backend

Já tentamos. Um servidor Python precisa de alguém que mantenha
dependências, renove certificado, acorde o serviço quando o free tier
hiberna e migre tudo quando o provedor muda as regras. Esse alguém é o
aluno, e o aluno vai defender e sair. Site estático continua no ar sozinho.

---

## O que o navegador faz

Ao abrir uma região, ele baixa o `.bin.gz`, descomprime, e monta o grafo de
fluxo: direção D8, acumulação por ordem topológica, e o grafo invertido em
formato CSR. Levou 39 ms para 77 mil células nos testes. A partir daí cada
clique custa de 0 a 2 ms, porque delimitar é só percorrer o grafo invertido
a partir do exutório, visitando apenas as células da própria bacia.

Depois da primeira visita o arquivo fica no cache, e o portal funciona sem
internet. Isso é útil de verdade para um técnico em campo.

---

## Limites conhecidos

- **D8 concentra o fluxo numa direção só.** Em encosta divergente, D-infinity
  daria contorno melhor. D8 foi escolhido porque deixa a delimitação exata e
  o grafo simples.
- **Terreno plano.** Em várzea e planície costeira o traçado sobre MDE de
  satélite é pouco confiável. O portal avisa quando a declividade do talvegue
  fica abaixo de 0,5 %, mas o certo é condicionar o MDE queimando a
  hidrografia oficial antes de processar.
- **Uma região por vez na memória.** Bacias que atravessam a borda do recorte
  saem truncadas. É por isso que o recorte tem que seguir divisor de águas.
- **Grupo hidrológico do solo é escolhido pelo usuário**, não lido de um mapa
  pedológico. É a maior fonte de incerteza que sobrou e o próximo módulo
  natural: cruzar com o mapa de solos e ponderar o CN também por ele.
- **A abstração inicial do SCS muda tudo.** Com λ = 0,20 e CN baixo, quase
  nada escoa e o SCS fica muito abaixo do racional. Com λ = 0,05 os dois
  convergem. Está exposto na interface justamente porque precisa de decisão
  consciente, e é bom tema de capítulo.

---

## Aviso

O portal entrega pré-dimensionamento. Não substitui projeto executivo nem
dispensa responsável técnico com ART.

---

## Publicando no GitHub Pages

Funciona, é gratuito e a estrutura deste pacote já está pronta para isso.
Três detalhes que quebram o site se passarem batido:

**A pasta precisa se chamar `docs/`.** O GitHub Pages só serve a partir do
root do repositório ou de uma pasta com esse nome exato. Por isso o site
está em `docs/` e não em `site/`.

**O arquivo `.nojekyll` precisa existir dentro de `docs/`.** Sem ele o
GitHub roda o Jekyll em cima da pasta, e o Jekyll ignora arquivos e pastas
que começam com underscore e mexe no que não devia. É um arquivo vazio, já
incluído. Ele não aparece em listagem normal de pasta, então ao arrastar
pelo navegador confira se subiu.

**O repositório precisa ser público.** Pages a partir de repositório privado
exige plano pago.

### Ativando

1. Crie o repositório e suba os arquivos preservando a estrutura de pastas.
   Pelo navegador, use "Add file" → "Upload files" e arraste a pasta inteira,
   não os arquivos soltos, senão tudo cai no root.
2. Settings → Pages → Source: "Deploy from a branch" → branch `main`,
   pasta `/docs` → Save.
3. Em um a dois minutos o endereço aparece como
   `https://SEU-USUARIO.github.io/NOME-DO-REPO/`.

O site usa caminhos relativos, então funciona normalmente nesse subcaminho.

### Quando o GitHub deixa de ser a melhor opção

O Git guarda todas as versões de todo arquivo, para sempre. Cada vez que
você regerar uma região de 3 MB, o repositório cresce mais 3 MB, mesmo
substituindo o arquivo. Com SC inteira e algumas rodadas de reprocessamento
isso passa fácil de 1 GB, que é o limite recomendado.

Enquanto forem poucas regiões, tudo bem. Quando virar o estado inteiro,
migre para Cloudflare Pages, que não versiona nada e não tem esse acúmulo.
O `index.html` não muda em nada nessa migração.
