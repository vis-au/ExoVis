// VIEW CONSTANTS

const MIN_PAA = 1;
const MAX_PAA = 14;
const MIN_SAX = 3;
const MAX_SAX = 14;

const MATRIX_CELL_SIZE = 45;
const MATRIX_WIDTH = (MAX_PAA) * (MATRIX_CELL_SIZE) * 2.5;
const MATRIX_HEIGHT = (MAX_SAX) * (MATRIX_CELL_SIZE);
const MATRIX_DUMMY_DOMAIN = [0.25, 89.11];

const PROGRESSION_BAR_SIZE = 100;
const PADDING_X = 30;
const PADDING_Y = 50;


// global variable that stores all cells selected by the user by either clicking a row, column
// or an individual cell.
const selectedCells = [];


// D3 COMPONENTS

// contains the root node of the visualization
const svg = d3.select("div#canvas").append("svg")
  .attr("width", MATRIX_WIDTH)
  .attr("height", MATRIX_HEIGHT)
  .style("background", "white");


// SCALES

// scale for columns (paa)
const scaleX = d3.scaleBand()
  .domain(d3.range(MIN_PAA, MAX_PAA + 1))
  .range([PADDING_X, MATRIX_WIDTH - PROGRESSION_BAR_SIZE - PADDING_X]);

// scale for rows (sax)
const scaleY = d3.scaleBand()
  .domain(d3.range(MIN_SAX, MAX_SAX + 1))
  .range([PADDING_Y + PROGRESSION_BAR_SIZE, MATRIX_HEIGHT - PADDING_Y]);

// scale for cell color
const color = d3.scaleSequentialLog(d3.interpolateYlGnBu)
  .domain(MATRIX_DUMMY_DOMAIN);

// scale for the column/row bars
const barScale = d3.scaleLinear().domain([0, 1]).range([0, PROGRESSION_BAR_SIZE]);

// size of each cell is xstep * ystep
const xStep = Math.abs(scaleX(MIN_PAA) - scaleX(MIN_PAA + 1));
const yStep = Math.abs(scaleY(MIN_SAX) - scaleY(MIN_SAX + 1));


// FUNCTIONS

function render(data_dictionary) {
  const data = getTransformedData(data_dictionary["matrix"]);
  const meanColumnProgress = data_dictionary["progression_column"].map(d => 1 - Math.abs(d));
  const meanRowProgress = data_dictionary["progression_row"].map(d => 1 - Math.abs(d));

  updateMatrix(data);
  renderBars(meanColumnProgress, meanRowProgress);
  renderAxes();
}

function notifyBackendSelectedCells() {
  const normalizedCells = selectedCells.map(cell => {
    return [cell.alpha, cell.omega];
  });

  window.eel.send_selected_cells(normalizedCells);
}

function toggleSelectedCell(cell) {
  const indexInSelected = getIndexInSelected(cell);

  if (indexInSelected === -1) {
    selectedCells.push(cell);
  } else {
    selectedCells.splice(indexInSelected, 1);
  }

  notifyBackendSelectedCells();
  updateSelectedCellStatus();
}

function selectColumn(columnIndex) {
  const omega = columnIndex + MIN_PAA;
  let allCellsAreSelected = true;
  let noCellsAreSelected = false;

  d3.range(MIN_SAX, MAX_SAX + 1).forEach(alpha => {
    const index = getIndexInSelected({ alpha, omega });
    allCellsAreSelected = allCellsAreSelected && index > -1;
    noCellsAreSelected = noCellsAreSelected && index === -1;
  });

  if (allCellsAreSelected || noCellsAreSelected) {
    d3.range(MIN_SAX, MAX_SAX + 1).forEach(alpha => {
      toggleSelectedCell({ alpha, omega });
    });
  } else {
    d3.range(MIN_SAX, MAX_SAX + 1).forEach(alpha => {
      const index = getIndexInSelected({ alpha, omega });
      if (index === -1) {
        toggleSelectedCell({ alpha, omega });
      }
    });
  }

  notifyBackendSelectedCells();
  updateSelectedCellStatus();
}

function selectRow(rowIndex) {
  const alpha = rowIndex + MIN_SAX;
  let allCellsAreSelected = true;
  let noCellsAreSelected = false;

  d3.range(MIN_PAA, MAX_PAA + 1).forEach(omega => {
    const index = getIndexInSelected({ alpha, omega });
    allCellsAreSelected = allCellsAreSelected && index > -1;
    noCellsAreSelected = noCellsAreSelected && index === -1;
  });

  if (allCellsAreSelected || noCellsAreSelected) {
    d3.range(MIN_PAA, MAX_PAA + 1).forEach(omega => {
      toggleSelectedCell({ alpha, omega });
    });
  } else {
    d3.range(MIN_PAA, MAX_PAA + 1).forEach(omega => {
      const index = getIndexInSelected({ alpha, omega });
      if (index === -1) {
        toggleSelectedCell({ alpha, omega });
      }
    });
  }

  notifyBackendSelectedCells();
  updateSelectedCellStatus();
}

function updateSelectedCellStatus() {
  svg.selectAll("rect.cell")
    .attr("stroke", d => getIndexInSelected(d) > -1 ? "white" : "none")
}

function getIndexInSelected(cell) {
  const { alpha, omega } = cell;

  return selectedCells.findIndex(selectedCell => {
    return selectedCell.alpha === alpha && selectedCell.omega === omega;
  });
}

function updateMatrix(data) {

  svg.selectAll("g.matrix").remove();

  const matrix = svg.append("g").attr("class", "matrix");

  const cell = matrix.selectAll("g.cell").data(data).join("g")
    .attr("class", "cell")
    .attr("transform", d => `translate(${scaleX(d.omega)}, ${scaleY(d.alpha)})`)
    .on("click", toggleSelectedCell)

  cell.append("title").text(d => d.error);

  cell.append("rect")
    .attr("class", "cell")
    .attr("width", xStep)
    .attr("height", yStep)
    .attr("fill", d => d.error === -1 ? "transparent" : color(d.error))
    .attr("stroke-width", 5);

  updateSelectedCellStatus();

  cell.append("text")
    .attr("class", "label")
    .attr("font-family", "sans-serif")
    .attr("font-size", 11)
    .attr("text-anchor", "middle")
    .attr("dx", xStep/2)
    .attr("dy", yStep * 0.66)
    .attr("fill", d => d.error < 9 ? (d.error === -1 ? "none" : "black" ): "white")
    .text(d => (d.error + "").slice(0, 4));
}


function renderAxes() {
  const xAxisGenerator = d3.axisBottom(scaleX);
  const yAxisGenerator = d3.axisLeft(scaleY);

  svg.selectAll("g.axis").remove();

  svg.append("g")
    .attr("class", "axis x")
    .attr("transform", `translate(0, ${MATRIX_HEIGHT - PADDING_Y})`)
    .call(xAxisGenerator);

  svg.append("g")
    .attr("class", "axis y")
    .attr("transform", `translate(${PADDING_X}, 0)`)
    .call(yAxisGenerator);
}

function renderBars(meanColumnProgress, meanRowProgress) {
  const barColor = "#afafaf";
  const barStroke = "#ccc";

  svg.selectAll("g.bars").remove();

  const barsX = svg.append("g").attr("class", "bars x");
  const barsY = svg.append("g").attr("class", "bars y");

  barsX.selectAll("rect.bar.x").data(meanColumnProgress).join("rect")
    .attr("class", "bar x")
    .attr("x", (d, i) => scaleX(i + 1))
    .attr("y", d => PROGRESSION_BAR_SIZE + PADDING_Y - barScale(d))
    .attr("width", xStep)
    .attr("height", d => barScale(d))
    .attr("fill", barColor)
    .attr("stroke", barStroke)
    .on("click", (d, i) => selectColumn(i));

  barsY.selectAll("rect.bar.y").data(meanRowProgress).join("rect")
    .attr("class", "bar y")
    .attr("x", scaleX(scaleX.domain()[MAX_PAA - MIN_PAA]) + xStep)
    .attr("y", (d, i) => scaleY(i + 3))
    .attr("width", d => barScale(d))
    .attr("height", yStep)
    .attr("fill", barColor)
    .attr("stroke", barStroke)
    .on("click", (d, i) => selectRow(i));

  renderIndicators();
}

function renderIndicators() {
  const barsX = svg.select("g.bars.x");
  const barsY = svg.select("g.bars.y");

  const indicators = [0.25, 0.5, 0.75, 1];

  const indicatorColumn = barsX.selectAll("g.indicator.column").data(indicators).join("g")
    .attr("class", "indicator column")
    .attr("transform", d => `translate(${PADDING_X},${PROGRESSION_BAR_SIZE + PADDING_Y - barScale(d)})`);

  indicatorColumn.append("line")
      .attr("x1", 0)
      .attr("x2", scaleX(MAX_PAA) + xStep - PADDING_X)
      .attr("y1", 0)
      .attr("y2", 0)
      .attr("stroke-dasharray", 2)
      .attr("stroke", "#333");

  indicatorColumn.append("text")
    .attr("dy", 4)
    .attr("dx", 4)
    .attr("x", scaleX(MAX_PAA) + xStep - PADDING_X)
    .attr("font-family", "sans-serif")
    .attr("font-size", 10)
    .text(d => parseInt(d * 100) + "%");


  const indicatorRow = barsY.selectAll("g.indicator.row").data(indicators).join("g")
    .attr("class", "indicator row")
    .attr("transform", d => `translate(${scaleX(scaleX.domain()[MAX_PAA - MIN_PAA]) + xStep + barScale(d)},${PROGRESSION_BAR_SIZE + PADDING_Y})`);

  indicatorRow.append("line")
      .attr("x1", 0)
      .attr("x2", 0)
      .attr("y1", 0)
      .attr("y2", scaleY(MAX_SAX) - scaleY(MIN_SAX) + yStep)
      .attr("stroke-dasharray", 2)
      .attr("stroke", "#333");

  indicatorRow.append("text")
    .attr("font-family", "sans-serif")
    .attr("font-size", 10)
    .attr("y", scaleY(MAX_SAX) - scaleY(MIN_SAX) + yStep + 10)
    .attr("text-anchor", "middle")
    .text(d => parseInt(d * 100));
}

/**
 * Simulate the stratified matrix data retrieved from the backend, including the progression on
 * columns and rows.
 */
function getDummyData() {
  const dummyData = [];
  const randomRowProgress = d3.range(MAX_SAX - MIN_SAX + 1).map(Math.random);
  const randomColumnProgress = d3.range(MAX_PAA - MIN_PAA + 1).map(Math.random);

  for (let paa = MIN_PAA; paa <= MAX_PAA; paa++) {
    for (let sax = MIN_SAX; sax <= MAX_SAX; sax++) {
      dummyData.push(Math.random() * MATRIX_DUMMY_DOMAIN[1] - MATRIX_DUMMY_DOMAIN[0]);
    }
  }

  const dummyObj = {
    matrix: dummyData,
    progression_column: randomColumnProgress,
    progression_row: randomRowProgress
  }

  return dummyObj;
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
          "omega": paa + 1,
          "error": arrayData[sax * MAX_SAX + paa],
          "progress": 100
        });
      }
    }
  }

  return transformedData;
}


// setInterval(() => {
  // render(getDummyData());
// }, 500);

window.eel.set_host("ws://localhost:8080");
window.eel.register_client("hello there");
window.eel.expose(render, "send_data");