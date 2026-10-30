-- Set pipe-table column widths from content length so long tables wrap instead of overflowing.
local function cell_len(cell)
  return utf8.len(pandoc.utils.stringify(cell.contents)) or 0
end

function Table(tbl)
  local n = #tbl.colspecs
  local maxlen, sumlen, cnt = {}, {}, {}
  for i = 1, n do maxlen[i], sumlen[i], cnt[i] = 1, 0, 0 end
  local function scan(rows)
    for _, row in ipairs(rows) do
      for i, cell in ipairs(row.cells) do
        if i <= n then
          local l = cell_len(cell)
          if l > maxlen[i] then maxlen[i] = l end
          sumlen[i] = sumlen[i] + l; cnt[i] = cnt[i] + 1
        end
      end
    end
  end
  scan(tbl.head.rows)
  for _, body in ipairs(tbl.bodies) do scan(body.body) end
  local w, total = {}, 0
  for i = 1, n do
    local avg = cnt[i] > 0 and sumlen[i] / cnt[i] or 1
    -- blend of mean and max length, with a floor so short columns stay legible
    w[i] = math.max(0.6 * avg + 0.4 * math.min(maxlen[i], 60), 6)
    total = total + w[i]
  end
  -- short tables (all cells short) keep natural widths
  local maxall = 0
  for i = 1, n do if maxlen[i] > maxall then maxall = maxlen[i] end end
  if total < 70 and maxall < 30 then return tbl end
  for i = 1, n do
    tbl.colspecs[i] = {tbl.colspecs[i][1], 0.98 * w[i] / total}
  end
  return tbl
end
