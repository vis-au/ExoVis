const MAX_PAA = 12;
const MAX_SAX = 14;

const MATRIX_DOMAIN = [0.25, 89.11];
const MATRIX_CELL_SIZE = 34;
const MATRIX_HISTOGRAM_HEIGHT = 50;
const MATRIX_HISTOGRAM_WIDTH = MATRIX_CELL_SIZE - 10;
const MATRIX_WIDTH = (MAX_PAA) * (MATRIX_CELL_SIZE) * 2.5;
const MATRIX_HEIGHT = (MAX_SAX) * (MATRIX_CELL_SIZE);
const MATRIX_HISTOGRAM_COLOR = "grey";

const svg = d3.select("div#canvas").append("svg")
  .attr("width", MATRIX_WIDTH)
  .attr("height", MATRIX_HEIGHT)
  .style("background", "white");

const barSize = 100;
const paddingX = 30;
const paddingY = 50;

const scaleX = d3.scaleLinear()
  .domain([0, MAX_PAA])
  .range([paddingX, MATRIX_WIDTH - barSize - paddingX]);

const scaleY = d3.scaleLinear()
  .domain([3, MAX_SAX + 1])
  .range([paddingY + barSize, MATRIX_HEIGHT - paddingY]);''

const xStep = scaleX(2) - scaleX(1);
const yStep = scaleY(2) - scaleY(1);

const color = d3.scaleSequentialLog(d3.interpolateYlGnBu)
  .domain(MATRIX_DOMAIN);


// let dummyData = getDummyData();
// sendRealData(dummyData);
// console.log(getTransformedData(dummyData))


function render(listData) {
  const data = getTransformedData(listData);
  const meanXProgress = d3.range(MAX_PAA).map(d => Math.random());
  const meanYProgress = d3.range(2, MAX_SAX).map(d => Math.random());

  renderAxes();
  updateMatrix(data);
  renderBars(meanXProgress, meanYProgress);
}

function updateMatrix(data) {

  svg.selectAll("g.matrix").remove();

  const matrix = svg.append("g").attr("class", "matrix");

  matrix.selectAll("rect.cell").data(data).join("rect")
    .attr("class", "cell")
    .attr("x", d => scaleX(d.omega))
    .attr("y", d => scaleY(d.alpha))
    .attr("width", xStep)
    .attr("height", yStep)
    .attr("fill", d => color(d.error))
    .attr("stroke", "none");

  matrix.selectAll("text.label").data(data).join("text")
    .attr("class", "label")
    .attr("font-family", "sans-serif")
    .attr("font-size", 11)
    .attr("text-anchor", "middle")
    .attr("x", d => scaleX(d.omega) + xStep/2)
    .attr("y", d => scaleY(d.alpha) + yStep * 0.66)
    .attr("fill", d => d.error < 30 ? "black" : "white")
    .text(d => (d.error + "").slice(0, 4));
}


function renderAxes() {
  const xAxisGenerator = d3.axisBottom(scaleX);
  const yAxisGenerator = d3.axisLeft(scaleY);

  svg.selectAll("g.axis").remove();

  const axisX = svg.append("g")
    .attr("class", "axis x")
    .attr("transform", `translate(0, ${MATRIX_HEIGHT - paddingY})`)
    .call(xAxisGenerator);

  const axisY = svg.append("g")
    .attr("class", "axis y")
    .attr("transform", `translate(${paddingX}, 0)`)
    .call(yAxisGenerator);

  axisX.append("text")
    .text("omega")
    .attr("x", paddingX + MATRIX_WIDTH/2)
    .attr("y", MATRIX_HEIGHT + paddingY);
}

const barScale = d3.scaleLinear().domain([0, 1]).range([0, barSize]);

function renderBars(meanXProgress, meanYProgress) {
  const barColor = "#afafaf";
  const barStroke = "#ccc";

  svg.selectAll("g.bars").remove();

  const barsX = svg.append("g").attr("class", "bars x");
  const barsY = svg.append("g").attr("class", "bars y");

  barsX.selectAll("rect.bar.x").data(meanXProgress).join("rect")
    .attr("class", "bar x")
    .attr("x", (d, i) => scaleX(i))
    .attr("y", d => barSize + paddingY - barScale(d))
    .attr("width", xStep)
    .attr("height", d => barScale(d))
    .attr("fill", barColor)
    .attr("stroke", barStroke);

  barsY.selectAll("rect.bar.y").data(meanYProgress).join("rect")
    .attr("class", "bar y")
    .attr("x", scaleX(scaleX.domain()[1]))
    .attr("y", (d, i) => scaleY(i + 3))
    .attr("width", d => barScale(d))
    .attr("height", yStep)
    .attr("fill", barColor)
    .attr("stroke", barStroke);

  barsX.append("text")
    .attr("transform")
}

render(getDummyData());


/**
 * Simulate the stratified matrix data retrieved from the backend
 */
function getDummyData() {
  const dummyData = [];

  for (let paa = 0; paa < MAX_PAA; paa++) {
    for (let sax = 0; sax < MAX_SAX; sax++) {
      dummyData.push(Math.random() * MATRIX_DOMAIN[1] - MATRIX_DOMAIN[0]);
    }
  }

  return dummyData;
}

/**
 * Transform stratified matrix data into tabular data that can be used by vega-lite. Assumes the
 * data to be stored row-wise in |paa| rows and |sax| columns.
 * @param {number[]} arrayData
 */
function getTransformedData(arrayData) {
  const transformedData = [];

  for (let paa = 0; paa < MAX_PAA; paa++) {
    for (let sax = 0; sax < MAX_SAX; sax++) {
      if (arrayData[sax * MAX_SAX + paa] !== undefined) {
        transformedData.push({
          "alpha": sax + 3,
          "omega": paa,
          "error": arrayData[sax * MAX_SAX + paa],
          "progress": 100
        });
      }
    }
  }

  return transformedData;
}


// data contains list of objects with 'cell' property for the parametrization and 'value' for the
// new entry
function sendDataToFrontend(next_chunk) {
  // transform cell into the format used in the vega-lite view
  const new_values = next_chunk.map(datum => {
    const paa = datum.cell[0];
    const sax = datum.cell[1];
    const value = datum.value;

    return { paa, sax, value };
  });

  // update matrix
  view.then(view => {
    view.insert("timeseries", new_values).run();
  });
}

window.eel.set_host("ws://localhost:8080");
window.eel.register_client("hello there");
// window.eel.expose(sendRealData, "send_data");
window.eel.expose(render, "send_data");