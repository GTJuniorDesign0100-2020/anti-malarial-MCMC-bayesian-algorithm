import React from 'react';
import { useState, useRef } from 'react';

/**
 * A generic collapsible accordion menu that can hold arbitrary content
 */
export default function Accordion({header, children}) {
  const [isCollapsed, setIsCollapsed] = useState(true);
  const hiddenContent = useRef(null);

  const contentStyle = {
    overflow: 'hidden',
    maxHeight: isCollapsed ? 0 : `${hiddenContent.current.scrollHeight}px`,
    transition: 'all 0.2s ease-out'
  };

  return (
    <div>
      <div onClick={() => setIsCollapsed(!isCollapsed)}>
        {header}
      </div>
      <div ref={hiddenContent} style={contentStyle}>
        {children}
      </div>
    </div>
  );
}
