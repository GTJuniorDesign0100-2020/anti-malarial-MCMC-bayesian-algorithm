import React from 'react';
import { useState, useMemo } from 'react';

/**
 * A generic table with sortable columns, based on this tutorial from Smashing:
 * https://www.smashingmagazine.com/2020/03/sortable-tables-react/
 *
 * columnNames: An array of names to use as the column names in the table
 * columnSortKeys: For each column, the key in each "item" to access when sorting
 * items: An array of items to display in the table
 * itemToTableRowFunc: A function that takes in 1 item in "items" and returns a
 * table row (including <tr>...</tr>)
 */
export default function SortableTable({
  columnNames,
  columnSortKeys,
  items,
  itemToTableRowFunc}) {
  const { sortedItems, requestSort, sortConfig } = useSortableData(items);
  const getClassNamesFor = (name) => {
    if (!sortConfig) {
      return;
    }
    return sortConfig.key === name ? sortConfig.direction : undefined;
  };

  console.log(sortedItems);
  console.log(sortedItems
    ? sortedItems.map(item => itemToTableRowFunc(item))
    : undefined);

  return (
    <table>
      <thead>
        <tr>
          {
            columnNames.map((columnName, i) => {
              return (
                <th key={columnName}>
                  <button
                    type="button"
                    onClick={() => requestSort(columnSortKeys[i])}
                    className={getClassNamesFor(columnSortKeys[i])}
                  >
                    {columnName}
                  </button>
                </th>
              );
            })
          }
        </tr>
      </thead>
      <tbody>
        {sortedItems
          ? sortedItems.map(item => itemToTableRowFunc(item))
          : undefined
        }
      </tbody>
    </table>
  );
};

const useSortableData = (items, config = null) => {
  const [sortConfig, setSortConfig] = useState(config);

  // Memoize to cache sort results, saving time
  const sortedItems = useMemo(() => {
    let sortableItems = [...items];
    if (sortConfig !== null) {
      sortableItems.sort((a, b) => {
        if (a[sortConfig.key] < b[sortConfig.key]) {
          return sortConfig.direction === 'ascending' ? -1 : 1;
        }
        if (a[sortConfig.key] > b[sortConfig.key]) {
          return sortConfig.direction === 'ascending' ? 1 : -1;
        }
        return 0;
      });
    }
    return sortableItems;
  }, [items, sortConfig]);

  const requestSort = key => {
    let direction = 'ascending';
    if (sortConfig && sortConfig.key === key && sortConfig.direction === 'ascending') {
      direction = 'descending';
    }
    setSortConfig({ key, direction });
  }

  return { sortedItems, requestSort, sortConfig };
};
