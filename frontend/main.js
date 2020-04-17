const MAX_PAA = 12;
const MAX_SAX = 14;

const MATRIX_DOMAIN = [0.25, 89.11];
const MATRIX_CELL_SIZE = 24;
const MATRIX_HISTOGRAM_HEIGHT = 50;
const MATRIX_HISTOGRAM_WIDTH = MATRIX_CELL_SIZE - 10;
const MATRIX_WIDTH = (MAX_PAA) * (MATRIX_CELL_SIZE);
const MATRIX_HEIGHT = (MAX_SAX) * (MATRIX_CELL_SIZE);
const MATRIX_HISTOGRAM_COLOR = "grey";


/**
 * Simulate the stratified matrix data retrieved from the backend
 */
function getDummyData() {
  const dummyData = [];

  for (let paa = 0; paa < MAX_PAA; paa++) {
    for (let sax = 0; sax < MAX_SAX; sax++) {
      dummyData.push(Math.random());
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
      transformedData.push({
        "alpha": sax,
        "omega": paa,
        "error": arrayData[sax * MAX_SAX + paa],
        "progress": 100
      });
    }
  }

  return transformedData;
}

// reference to the vega embed view, is set by the vegaEmbed promise on rendering.
let view = null;

const paaScale = {
  "domain": [0, MAX_PAA]
};
const saxScale = {
  "domain": [MAX_SAX, 0]
};

function sendRealData(realData) {
  const spec = {
    "$schema": "https://vega.github.io/schema/vega-lite/v4.json",
    "data": {
      values: getTransformedData(realData)
    },
    "vconcat": [
      {
        "mark": "bar",
        "width": MATRIX_WIDTH,
        "height": MATRIX_HISTOGRAM_HEIGHT,
        "encoding": {
          "x": {
            "bin": {
              "maxbins": MAX_PAA
            },
            "field": "omega",
            "type": "quantitative",
            "axis": null,
            scale: paaScale
          },
          "y": {
            "aggregate": "mean",
            "field": "progress",
            "type": "quantitative",
            "scale": {
              "domain": [0, 100]
            },
            "axis": null,
          },
          "color": {
            "value": MATRIX_HISTOGRAM_COLOR
          },
          "size": {
            "value": MATRIX_HISTOGRAM_WIDTH
          }
        }
      },
      {
        "hconcat": [
          {
            "layer": [
              {
                "mark": "square",
                "width": MATRIX_WIDTH,
                "height": MATRIX_HEIGHT,
                "encoding": {
                  "x": {
                    "field": "omega",
                    "type": "quantitative",
                    scale: paaScale,
                    "axis": null,
                  },
                  "y": {
                    "field": "alpha",
                    "type": "quantitative",
                    scale: saxScale,
                    "axis": null,
                  },
                  "size": {
                    "value": 600
                  },
                  "color": {
                    "field": "error",
                    "type": "quantitative",
                    "scale": {
                      "scheme": "yellowgreenblue",
                      domain: MATRIX_DOMAIN,
                      "type": "log",
                      "base": 10
                    },
                    "legend": null
                  },
                  "opacity": {
                    "value": 1
                  }
                }
              },
              // {
              //   "mark": "text",
              //   "width": mainWidth,
              //   "height": mainHeight,
              //   "encoding": {
              //     "x": {
              //       "field": "omega",
              //       "type": "quantitative",
              //       scale: paaScale
              //     },
              //     "y": {
              //       "field": "alpha",
              //       "type": "quantitative",
              //       scale: saxScale
              //     },
              //     "text": {
              //       "field": "error",
              //       "format": "2f"
              //     }
              //   }
              // }
            ]
          },
          {
            "mark": "bar",
            "height": MATRIX_HEIGHT,
            "width": MATRIX_HISTOGRAM_HEIGHT,
            "encoding": {
              "y": {
                "bin": {
                  "maxbins": MAX_SAX
                },
                "field": "alpha",
                "type": "quantitative",
                "axis": null,
                scale: saxScale
              },
              "x": {
                "aggregate": "mean",
                "field": "progress",
                "type": "quantitative",
                "scale": {
                  "domain": [0, 100]
                },
                "axis": null,
              },
              "size": {
                "value": MATRIX_HISTOGRAM_WIDTH
              },
              "color": {
                "value": MATRIX_HISTOGRAM_COLOR
              },
            }
          }
        ]
      }
    ],
    "config": {
      "concat": {
        "spacing": 0
      },
      "view": {
        "stroke": "transparent"
      },
      "padding": 0
    }
  };

  view = vegaEmbed("div#canvas", spec)
    .then(res => res.view);
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
window.eel.expose(sendRealData, "send_data");