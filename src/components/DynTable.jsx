import React from 'react';
import JSZip from 'jszip';
import {saveAs} from 'file-saver';
import SortableTable from './SortableTable';

export default function DynTable(props) {
  function renderResultRow(resultData) {
    const {date, inputFilename, status, results} = resultData;
    const csvFileText = results.output_file_text;
    return (
      <tr key={date.toISOString()}>
        <td>{date.toLocaleDateString()}</td>
        <td>{date.toLocaleTimeString()}</td>
        <td>{inputFilename}</td>
        <td>{status}</td>
        <td>{getResultDownloadLinks(csvFileText)}</td>
      </tr>
    );
  }

  return (
    <div>
      <h1 align='center'>Results</h1>

      <SortableTable
        columnNames={['Date', 'Time', 'Input File Name', 'Status', 'Output']}
        columnSortKeys={['date', 'date', 'inputFilename', 'status', '']}
        items={
          // Reverse so newest results appear at the top
          Object.values(props.data).reverse()
        }
        itemToTableRowFunc={renderResultRow}
      />
    </div>
  );
}

function getResultDownloadLinks(csvFileText) {
  if (!csvFileText) {
    return '';
  }

  const probability_filename = 'probability_of_recrudescence.csv';

  let extra_filenames = Object.keys(csvFileText);
  // Remove probability file from those included in .zip
  extra_filenames.splice(extra_filenames.indexOf(probability_filename), 1);

  return (
    <div>
      <CSVDownloadLink
        key={probability_filename}
        csvFileName={probability_filename}
        csvFileText={csvFileText[probability_filename]}
      />
      <ZipDownloadLink
        fileNames={extra_filenames}
        fileContentDict={csvFileText}
        downloadName="advanced_stats"
      />
    </div>
  );
}

const ZipDownloadLink = ({fileNames, fileContentDict, downloadName}) => {
  let zip = new JSZip();
  for (let filename of fileNames) {
    const csvText = fileContentDict[filename];
    zip.file(filename, csvText);
  }

  return (
    <p><a
      href="#"
      onClick={(evt) => {
        evt.preventDefault();
        zip.generateAsync({type: 'blob'}).then(
          zipContent => saveAs(zipContent, `${downloadName}.zip`)
        );
      }}>
        {downloadName}.zip
    </a></p>
  );
};

const CSVDownloadLink = ({csvFileName, csvFileText}) => {
  const csvText = `data:text/csv;charset=utf-8,${csvFileText}`;
  const csvUri = encodeURI(csvText);
  return (
    <p><a href={csvUri} download={csvFileName}>{csvFileName}</a></p>
  );
};
